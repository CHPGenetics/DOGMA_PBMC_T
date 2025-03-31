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

meta <- meta[meta$condition %in% c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'),]
meta$label <- factor(meta$condition, levels = c('Act_IL1B_IL23_PGE2', 'Act_IL1B_IL23'))
meta$cell_type <-factor(meta$celltype_updated, levels = levels)
meta$replicate <- meta$sample
table(meta$cell_type)
rna <- rna[,rownames(meta)]
adt <- adt[,rownames(meta)]
peak <- peak[,rownames(meta)]

rna_mat = to_pseudobulk(rna, meta = meta)
peak_mat = to_pseudobulk(peak, meta = meta)

eQTA = function(celltype){
  print(celltype)
  rna = varianceStabilizingTransformation(as.matrix(rna_mat[[celltype]]))
  atac = varianceStabilizingTransformation(as.matrix(peak_mat[[celltype]]))
  
  saveRDS(rna, file = paste0('../output/eQTA_updated/matrix_2/RNA/', celltype, '.RDS'))
  saveRDS(atac, file = paste0('../output/eQTA_updated/matrix_2/ATAC/', celltype, '.RDS'))
}

for (i in 1:20){
  eQTA(levels[i])
}