rm(list=ls())
gc()

# Load Packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)

# Parameters
DATA_PATH <- "/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/01_Th1_Th17/"

# Load Data
tcell_dogma <- readRDS(paste0(DATA_PATH, "tcell_annotated_updated.RDS"))

# DE
Th_dogma <- subset(tcell_dogma, subset = celltype_updated %in% c("", ""))
DefaultAssay(Th_dogma) <- 'RNA'



# DA
DefaultAssay(Th_dogma) <- 'peaks'
Idents(Th_dogma) <- "celltype_updated"
da_peaks <- FindMarkers(
  object = Th_dogma,
  ident.1 = c("CD4+ Memory (Resting) - Th1", "CD4+ Memory (Activated) - Th1"), 
  ident.2 = c("CD4+ Memory (Resting) - Th17", "CD4+ Memory (Activated) - Th17"),
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)







