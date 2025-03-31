library(Libra)
library(data.table)
library(parallel)
library(dplyr)
library(magrittr)
library(tibble)
library(forcats)
library(DESeq2)
library(edgeR)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

meta <- read.csv('../output/tcell_annotated_updated.csv', row.names = 'X')
rna <- readRDS('../output/pseudo_bulk/rna.RDS')
adt <- readRDS('../output/pseudo_bulk/adt.RDS')
peak <- readRDS('../output/pseudo_bulk/peak.RDS')
peak_annotation <- readRDS('../output/pseudo_bulk/peak_annotation.RDS')

meta <- meta[meta$condition %in% c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'),]
meta$label <- factor(meta$condition, levels = c('Act_IL1B_IL23_PGE2', 'Act_IL1B_IL23'))
meta$cell_type <- meta$celltype_updated
meta$replicate <- meta$sample
meta$cell_type <- gsub(' (Activated)', '', meta$cell_type, fixed = T)
meta$cell_type <- gsub(' (Resting)', '', meta$cell_type, fixed = T)
table(meta$cell_type)
rna <- rna[,rownames(meta)]
adt <- adt[,rownames(meta)]
peak <- peak[,rownames(meta)]

pseudobulk_de = function(input = input, 
                         meta = meta, 
                         de_family = de_family,
                         de_method = de_method,
                         de_type = de_type) {
  # check args
  if (de_method == 'limma') {
    if (de_type != 'voom') {
      # change default type to use
      de_type = 'trend'  
    }
  }
  
  # get the pseudobulks list
  pseudobulks = to_pseudobulk(
    input = input,
    meta = meta
  )
  
  results = purrr::map(pseudobulks, function(x) {
    # create targets matrix
    targets = data.frame(group_sample = colnames(x)) %>%
      mutate(group = gsub(".*\\:", "", group_sample), sample = gsub("\\:.*", "", group_sample))
    ## optionally, carry over factor levels from entire dataset
    if (is.factor(meta$label)) {
      targets$group %<>% factor(levels = levels(meta$label))
    }
    if (n_distinct(targets$group) > 2)
      return(NULL)
    # create design
    design = model.matrix(~ group + sample, data = targets)
    
    DE = switch(de_method,
                edgeR = {
                  tryCatch({
                    y = DGEList(counts = x, group = targets$group) %>%
                      calcNormFactors(method = 'TMM') %>%
                      estimateDisp(design)
                    test = switch(de_type,
                                  QLF = {
                                    fit = glmQLFit(y, design)
                                    test = glmQLFTest(fit, coef = -1)
                                  },
                                  LRT = {
                                    fit = glmFit(y, design = design)
                                    test = glmLRT(fit)
                                  })
                    res = topTags(test, n = Inf) %>%
                      as.data.frame() %>%
                      rownames_to_column('gene') %>%
                      # flag metrics in results
                      mutate(de_family = 'pseudobulk',
                             de_method = de_method,
                             de_type = de_type)
                  }, error = function(e) {
                    message(e)
                    data.frame()
                  })
                },
                DESeq2 = {
                  tryCatch({
                    dds = DESeqDataSetFromMatrix(countData = x,
                                                 colData = targets,
                                                 design = ~ sample + group)
                    dds = switch(de_type,
                                 Wald = {
                                   dds = try(DESeq(dds,
                                                   test = 'Wald',
                                                   fitType = 'parametric',
                                                   sfType = 'poscounts',
                                                   betaPrior = F))
                                 },
                                 LRT = {
                                   dds = try(DESeq(dds,
                                                   test = 'LRT',
                                                   reduced = ~ sample,
                                                   fitType = 'parametric',
                                                   sfType = 'poscounts',
                                                   betaPrior = F))
                                 }
                    )
                    res = results(dds)
                    # write
                    res = as.data.frame(res) %>%
                      mutate(gene = rownames(x)) %>%
                      # flag metrics in results
                      mutate(de_family = 'pseudobulk',
                             de_method = de_method,
                             de_type = de_type)
                  }, error = function(e) {
                    message(e)
                    data.frame()
                  })
                },
                limma = {
                  tryCatch({
                    x = switch(de_type,
                               trend = {
                                 trend_bool = T
                                 dge = DGEList(as.matrix(x), group = targets$group)
                                 dge = calcNormFactors(dge)
                                 x = new("EList")
                                 x$E = cpm(dge, log = TRUE, prior.count = 3)
                                 x
                               },
                               voom = {
                                 counts = all(as.matrix(x) %% 1 == 0)
                                 if (counts) {
                                   trend_bool = F
                                   x = voom(as.matrix(x), design)
                                   x
                                 }
                               })
                    # get fit
                    fit = lmFit(x, design) %>%
                      eBayes(trend = trend_bool, robust = trend_bool)
                    # format the results
                    res = fit %>%
                      # extract all coefs except intercept
                      topTable(number = Inf, coef = -1) %>%
                      rownames_to_column('gene') %>%
                      # flag metrics in results
                      mutate(
                        de_family = 'pseudobulk',
                        de_method = de_method,
                        de_type = de_type)
                  }, error = function(e) {
                    message(e)
                    data.frame()
                  })
                }
    )
  })
  results %<>% bind_rows(.id = 'cell_type')
}

run_de_paired <- function(input, 
                          meta = NULL, 
                          de_family = 'pseudobulk',
                          de_method = 'edgeR',
                          de_type = 'LRT') {
  DE <- pseudobulk_de(input = input, 
                      meta = meta, 
                      de_family = de_family,
                      de_method = de_method,
                      de_type = de_type)
  # clean up the output
  suppressWarnings(
    colnames(DE) %<>%
      fct_recode('p_val' = 'p.value',  ## DESeq2
                 'p_val' = 'pvalue',  ## DESeq2
                 'p_val' = 'p.value',  ## t/wilcox
                 'p_val' = 'P.Value',  ## limma
                 'p_val' = 'PValue'  , ## edgeR
                 'p_val_adj' = 'padj', ## DESeq2/t/wilcox
                 'p_val_adj' = 'adj.P.Val',      ## limma
                 'p_val_adj' = 'FDR',            ## edgeER
                 'avg_logFC' = 'log2FoldChange', ## DESEeq2
                 'avg_logFC' = 'logFC', ## limma/edgeR
                 'avg_logFC' = 'avg_log2FC' # Seurat V4
      )
  ) %>%
    as.character()
  
  DE %<>%
    # calculate adjusted p values
    group_by(cell_type) %>%
    mutate(p_val_adj = p.adjust(p_val, method = 'BH')) %>%
    # make sure gene is a character not a factor
    mutate(gene = as.character(gene)) %>%
    # invert logFC to match Seurat level coding
    mutate(avg_logFC = avg_logFC * -1) %>%
    dplyr::select(cell_type,
                  gene,
                  avg_logFC,
                  p_val,
                  p_val_adj,
                  de_family,
                  de_method,
                  de_type
    ) %>%
    ungroup() %>%
    arrange(cell_type, gene)
}

DE_rna <- run_de_paired(
  rna,
  meta = meta,
  de_family = "pseudobulk",
  de_method = "DESeq2",
  de_type = "LRT"
)

fwrite(DE_rna, file = '../output/pseudo_bulk/DE_rna_DESeq2_paired_merged.txt', col.names = T, row.names = F, quote = F)

DE_adt <- run_de_paired(
  adt,
  meta = meta,
  de_family = "pseudobulk",
  de_method = "DESeq2",
  de_type = "LRT"
)

fwrite(DE_adt, file = '../output/pseudo_bulk/DE_adt_DESeq2_paired_merged.txt', col.names = T, row.names = F, quote = F)

DE_peak <- run_de_paired(
  peak,
  meta = meta,
  de_family = "pseudobulk",
  de_method = "DESeq2",
  de_type = "LRT"
)

DE_peak <- cbind(DE_peak, peak_annotation[DE_peak$gene,c(2,4,5,6,8)])

fwrite(DE_peak, file = '../output/pseudo_bulk/DE_peak_DESeq2_paired_merged.txt', col.names = T, row.names = F, quote = F)

