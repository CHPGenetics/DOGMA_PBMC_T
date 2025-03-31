library(Libra)
library(data.table)
library(parallel)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

meta <- read.csv('../output/tcell_annotated_updated.csv', row.names = 'X')
rna <- readRDS('../output/pseudo_bulk/rna.RDS')
adt <- readRDS('../output/pseudo_bulk/adt.RDS')
peak <- readRDS('../output/pseudo_bulk/peak.RDS')
peak_annotation <- readRDS('../output/pseudo_bulk/peak_annotation.RDS')

meta <- meta[meta$condition %in% c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'),]
meta$condition <- factor(meta$condition, levels = c('Act_IL1B_IL23_PGE2', 'Act_IL1B_IL23'))
rna <- rna[,rownames(meta)]
adt <- adt[,rownames(meta)]
peak <- peak[,rownames(meta)]

DE_rna <- run_de(
  rna,
  meta = meta,
  replicate_col = "sample",
  cell_type_col = "celltype_updated",
  label_col = "condition",
  min_cells = 3,
  min_reps = 2,
  min_features = 0,
  de_family = "pseudobulk",
  de_method = "DESeq2",
  de_type = "LRT",
  n_threads = 48
)

fwrite(DE_rna, file = '../output/pseudo_bulk/DE_rna_DESeq2.txt', col.names = T, row.names = F, quote = F)

DE_adt <- run_de(
  adt,
  meta = meta,
  replicate_col = "sample",
  cell_type_col = "celltype_updated",
  label_col = "condition",
  min_cells = 3,
  min_reps = 2,
  min_features = 0,
  de_family = "pseudobulk",
  de_method = "DESeq2",
  de_type = "LRT",
  n_threads = 48
)

fwrite(DE_adt, file = '../output/pseudo_bulk/DE_adt_DESeq2.txt', col.names = T, row.names = F, quote = F)

DE_peak <- run_de(
  peak,
  meta = meta,
  replicate_col = "sample",
  cell_type_col = "celltype_updated",
  label_col = "condition",
  min_cells = 3,
  min_reps = 2,
  min_features = 0,
  de_family = "pseudobulk",
  de_method = "DESeq2",
  de_type = "LRT",
  n_threads = 48
)
DE_peak <- cbind(DE_peak, peak_annotation[DE_peak$gene,c(2,4,5,6,8)])

fwrite(DE_peak, file = '../output/pseudo_bulk/DE_peak_DESeq2.txt', col.names = T, row.names = F, quote = F)

