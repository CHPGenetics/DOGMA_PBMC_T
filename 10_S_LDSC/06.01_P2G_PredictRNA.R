rm(list=ls())
gc()

# Load Packages
library(Seurat)
library(Signac)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(parallel)
library(readxl)

# Parameters
DATA_PATH <- "/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/04_Global_Analysis/04_P2G_Shuffle/"

tcell_dogma <- readRDS(paste0(DATA_PATH, "tcell_annotated_updated.RDS"))
tcell_dogma <- subset(tcell_dogma, subset = condition %in% c("Act_IL1B_IL23"))

# Select Genes
top_gene <- read_xlsx("/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/04_Global_Analysis/03_P2G/Colocalization.xlsx", 1)
gene_list <- strsplit(top_gene$gene, "\\|")
gene_vector <- unique(unlist(gene_list))

rna.mat <- tcell_dogma@assays$RNA@data[gene_vector,]
saveRDS(rna.mat, paste0(WORK_PATH, "rna.rds"))

set.seed(20240404)
shuffled_matrix <- apply(rna.mat, 2, function(col) sample(col, length(col)))
shuffled_matrix <- matrix(shuffled_matrix, nrow = nrow(rna.mat), 
                          dimnames = dimnames(rna.mat))
saveRDS(shuffled_matrix, paste0(WORK_PATH, "rna.shuffle.rds"))

# Peak
peak.mat <- tcell_dogma@assays$peaks@data
saveRDS(peak.mat, paste0(WORK_PATH, "peak.rds"))

meta <- tcell_dogma@meta.data
saveRDS(meta, paste0(WORK_PATH, "meta.rds"))

tcell_dogma@assays <- tcell_dogma@assays[c("peaks")]
gc()

# quantify gene activity
DefaultAssay(tcell_dogma) <- "peaks"
gene.activities <- GeneActivity(tcell_dogma, features = gene_vector)
saveRDS(gene.activities, paste0(WORK_PATH, "rna.predicted.rds"))




