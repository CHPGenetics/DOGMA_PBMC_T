rm(list=ls())
gc()

# Load Packages
library(optparse)
library(data.table)
library(readxl)
library(tidyverse)
library(Seurat)
library(Signac)

## Input parameters
args_list = list(
  make_option("--celltype", type = "character", default = NULL,
              help = "INPUT: celltype to compare",
              metavar = "character"))

opt_parser = OptionParser(option_list = args_list)
opt = parse_args(opt_parser)

# Path Parameters
DATA_PATH <- "/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/pseudo_bulk_updated/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Input/Pseudo_Celltype_DA/"

# Load Data
tcell_dogma <- readRDS("/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/tcell_annotated_updated.RDS")
all_clusters <- unique(tcell_dogma$celltype_updated)

Idents(tcell_dogma) <- "condition"
tcell_dogma <- subset(x = tcell_dogma, idents = c("Act_IL1B_IL23_PGE2"))

# DA
Idents(tcell_dogma) <- "celltype_updated"
DefaultAssay(tcell_dogma) <- 'peaks'
DA_Markers <- FindMarkers(tcell_dogma, ident.1 = opt$celltype, ident.2 = NULL)
DA_Markers$Peak <- rownames(DA_Markers)
write.table(DA_Markers,
            file = paste0(WORK_PATH, "DA_", gsub(" ", "_", opt$celltype), ".txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)

