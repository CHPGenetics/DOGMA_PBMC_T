library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(data.table)
library(Matrix)
library(dplyr)
library(stringr)
library(future)
library(harmony)
library(parallel)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')

load('../data/data_fragments.RData')

ref_peaks <- readRDS('~/RWorkSpace/DOGMA-seq/DOGMA_analysis/output/peak_calling/ref_peaks.RDS')

plan('multicore')
options(future.globals.maxSize= 60*1024^3)

macs2_counts <- FeatureMatrix(
  fragments = fragments,
  features = ref_peaks,
  cells = colnames(data)
)

save(macs2_counts, file = '../data/macs2_counts_ref.RData')