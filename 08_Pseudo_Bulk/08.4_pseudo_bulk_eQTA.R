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
library(EnsDb.Hsapiens.v86)

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
meta = meta[meta$cell_type %in% c('CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
           'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1',
           'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17',
           'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh',
           'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other'),]
meta$cell_type <- "Pseudo-bulk"
table(meta$cell_type)
rna <- rna[,rownames(meta)]
adt <- adt[,rownames(meta)]
peak <- peak[,rownames(meta)]

rna_mat = to_pseudobulk(rna, meta = meta)
peak_mat = to_pseudobulk(peak, meta = meta)

rna = varianceStabilizingTransformation(as.matrix(rna_mat[[1]]))
atac = varianceStabilizingTransformation(as.matrix(peak_mat[[1]]))

table(colnames(atac) == colnames(rna))

gene <- as.data.frame(genes(EnsDb.Hsapiens.v86))
gene_10 <- fread('~/RWorkSpace/DOGMA-seq/Duerr_20220629_DOGMAseq-4_arc/filtered_feature_bc_matrix/features.tsv.gz')
gene_10 <- gene_10[gene_10$V3 == 'Gene Expression',]

table(rownames(rna) %in% gene$gene_name)#22586
table(rownames(rna) %in% gene_10$V2)#36591

colnames(gene_10) <- c('gene_id', 'gene_name', 'type', 'chr', 'start', 'end')

rna_loc <- gene_10

table(rownames(rna) %in% rna_loc$gene_name)#36591

shared_gene <- intersect(rownames(rna), rna_loc$gene_name)
rna_sub <- rna[shared_gene,]
rna_loc <- rna_loc[!duplicated(rna_loc$gene_name),]
rna_loc <- tibble::column_to_rownames(rna_loc, 'gene_name')
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
save(me, file = '../plots/eQTA/pseudobulk_eQTA.RData')