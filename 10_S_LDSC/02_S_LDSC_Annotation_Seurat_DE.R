rm(list=ls())
gc()

# Load Packages
library(data.table)
library(readxl)
library(tidyverse)

# Parameters
DATA_PATH <- "/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/pseudo_bulk_updated/"
REF_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Input/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Input/Seurat_DE/"

# Load Data
# tcell_dogma <- readRDS("/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/tcell_annotated_updated.RDS")
all_clusters <- unique(tcell_dogma$celltype_updated)

# ldcts
LDSC_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Output/"
cluster_file <- gsub(" ", "_", all_clusters)
path_list <- paste0(LDSC_PATH, "Seurat_DE/", cluster_file, "_", ",", LDSC_PATH, "Annotation_Gene/ref_")
ldcts <- data.frame(cluster_file, path_list)
write.table(ldcts, row.names = F, col.names = F, quote = F, sep = "\t",
            file = paste0(LDSC_PATH, "Seurat_DE/Result/ldcts"))
