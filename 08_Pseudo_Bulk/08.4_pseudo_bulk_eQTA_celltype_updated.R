library(Libra)
library(data.table)
library(parallel)
library(dplyr)
library(magrittr)
library(tibble)
library(forcats)
library(DESeq2)
library(edgeR)
library(MatrixEQTL)
library(stringr)
library(ArchR)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

meta <- read.csv('../output/tcell_annotated_updated.csv', row.names = 'X')
rna <- readRDS('../output/pseudo_bulk/rna.RDS')
adt <- readRDS('../output/pseudo_bulk/adt.RDS')
peak <- readRDS('../output/pseudo_bulk/peak.RDS')
peak_annotation <- readRDS('../output/pseudo_bulk/peak_annotation.RDS')

levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
           'CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
           'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1',
           'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17',
           'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh',
           'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other',
           'CD8+ Naive (Resting)', 'CD8+ Naive (Activated)',
           'CD8+ Regulatory',
           'CD8+ Memory (Resting)', 'CD8+ Memory (Activated)',
           'MAITs (Resting)', 'MAITs (Activated)', 'Gamma Delta'
)

reform_p2g <- function(p2g){
  peak <- metadata(p2g)$peakSet
  gene <- metadata(p2g)$geneSet
  p2g_df <- as.data.frame(p2g)
  
  peak_sub <- peak[p2g_df$idxATAC]
  gene_sub <- gene[p2g_df$idxRNA]
  peak_df <- data.frame(chr = seqnames(peak_sub), start = start(peak_sub), end = end(peak_sub))
  gene_df <- data.frame(chr = seqnames(gene_sub), pos = start(gene_sub), gene = gene_sub$name)
  peak_df$chr_pos <- paste(peak_df$chr, peak_df$start, peak_df$end, sep = '-')
  
  p2g_df$peak_chr <- peak_df$chr
  p2g_df$peak_start <- peak_df$start
  p2g_df$peak_end <- peak_df$end
  p2g_df$peak_chr_pos <- peak_df$chr_pos
  p2g_df$gene_chr <- gene_df$chr
  p2g_df$gene_pos <- gene_df$pos
  p2g_df$gene <- gene_df$gene
  length(unique(p2g_df$peak_chr_pos))#188999
  length(unique(p2g_df$gene))#25736
  p2g_df$distance <- pmin(abs(p2g_df$peak_start - p2g_df$gene_pos), abs(p2g_df$peak_end - p2g_df$gene_pos))
  p2g_df$distance[p2g_df$peak_start <= p2g_df$gene_pos & p2g_df$gene_pos <= p2g_df$peak_end] <- 0
  
  #p2g_df <- p2g_df[,-c(1:2)]
  return(p2g_df)
}

p2g_celltype_link <- list()
for(i in 1:9){
  celltype <- levels[i]
  name <- paste0('0',i,'_',celltype)
  p2g <- readRDS(paste0('../output/ArchR/DOGMA_filtered_multiome/', name, '/p2g_fdr_0.05.RDS'))
  p2g_celltype_link[[celltype]] <- reform_p2g(p2g)
}
for(i in 10:20){
  celltype <- levels[i]
  name <- paste0(i,'_',celltype)
  p2g <- readRDS(paste0('../output/ArchR/DOGMA_filtered_multiome/', name, '/p2g_fdr_0.05.RDS'))
  p2g_celltype_link[[celltype]] <- reform_p2g(p2g)
}

rna_loc_list <- lapply(p2g_celltype_link, function(p2g){
  loc <- p2g[,c('gene', 'gene_chr', 'gene_pos')]
  loc <- loc[!duplicated(loc$gene),]
  return(loc)
  })

rna_loc <- Reduce(rbind, rna_loc_list)
rna_loc <- rna_loc[!duplicated(rna_loc$gene),]
rownames(rna_loc) <- rna_loc$gene
colnames(rna_loc) <- c('gene_name', 'chr', 'start')

meta$label <- factor(meta$condition, levels = c('Act_IL1B_IL23_PGE2', 'Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2_TGFB', 'Act_IL1B_IL23_TGFB'))
meta$cell_type <- factor(meta$celltype_updated, levels = levels)
meta$replicate <- meta$sample
table(meta$cell_type)
rna <- rna[,rownames(meta)]
adt <- adt[,rownames(meta)]
peak <- peak[,rownames(meta)]

rna_mat = to_pseudobulk(rna, meta = meta)
peak_mat = to_pseudobulk(peak, meta = meta)

eQTA = function(celltype, rna_loc){
  print(celltype)
  rna = varianceStabilizingTransformation(as.matrix(rna_mat[[celltype]]))
  atac = varianceStabilizingTransformation(as.matrix(peak_mat[[celltype]]))
  
  table(colnames(atac) %in% colnames(rna))
  
  shared <- intersect(colnames(atac), colnames(rna))
  
  atac <- atac[,shared]
  rna <- rna[,shared]
  table(colnames(atac) == colnames(rna))
  
  table(rownames(rna) %in% rna_loc$gene_name)#36591
  
  shared_gene <- intersect(rownames(rna), rna_loc$gene_name)
  rna_sub <- rna[shared_gene,]
  rna_loc <- rna_loc[shared_gene,]
  table(rownames(rna_sub) == rownames(rna_loc))
  
  atac_loc <- as.data.frame(str_split_fixed(rownames(atac), '-', 3))
  colnames(atac_loc) = c('CHR', 'START', 'END')
  rownames(atac_loc) <- rownames(atac)
  
  rnapos <- data.frame(geneid = rownames(rna_loc), chr = rna_loc$chr, left = rna_loc$start, right = rna_loc$start)
  atacpos <- data.frame(snpid = rownames(atac_loc), chr = atac_loc$CHR, pos = round((as.numeric(atac_loc$START) + as.numeric(atac_loc$END))/2))
  
  #######
  atac_set = SlicedData$new()
  atac_set$CreateFromMatrix(atac)
  
  rna_set =  SlicedData$new()
  rna_set$CreateFromMatrix(rna_sub)
  
  saveRDS(rna_sub, file = paste0('../output/eQTA_updated/matrix/RNA/', celltype, '.RDS'))
  saveRDS(atac, file = paste0('../output/eQTA_updated/matrix/ATAC/', celltype, '.RDS'))
  
  library(parallel)
  me = Matrix_eQTL_main(
    snps = atac_set, 
    gene = rna_set, 
    cvrt = SlicedData$new(),
    pvOutputThreshold = 0,
    useModel = modelLINEAR, 
    errorCovariance = numeric(), 
    verbose = TRUE,
    output_file_name.cis = "",
    pvOutputThreshold.cis = 0.05,
    snpspos = atacpos, 
    genepos = rnapos,
    cisDist = 1e6,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
  saveRDS(me, file = paste0('../output/eQTA_updated/pseudobulk_eQTA_', celltype, '.RDS'))
}

for (i in 1:20){
  eQTA(levels[i], rna_loc)
}