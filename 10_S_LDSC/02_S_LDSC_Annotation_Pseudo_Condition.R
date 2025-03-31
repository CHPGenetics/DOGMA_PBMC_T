rm(list=ls())
gc()

# Load Packages
library(data.table)
library(readxl)
library(tidyverse)

# Parameters
DATA_PATH <- "/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/pseudo_bulk_updated/"
REF_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Input/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Input/Pseudo_Condition/"

# Load Data
# tcell_dogma <- readRDS("/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/tcell_annotated_updated.RDS")
all_clusters <- unique(tcell_dogma$celltype_updated)
DE <- read_xlsx(paste0(DATA_PATH, "DE_rna_DESeq2_paired_new.xlsx"), 1)
DA <- read_xlsx(paste0(DATA_PATH, "DE_peak_DESeq2_paired_new.xlsx"), 1)

# Load Reference
Genes <- read.table(paste0(REF_PATH, "Gene_coord_all.txt"), header = T)
Peaks <- read.table(paste0(REF_PATH, "Peak2_coord_all.txt"), header = T)
DE <- DE %>% filter(gene %in% Genes$GENE)
DA <- DA %>% filter(peak %in% Peaks$GENE)

# DE
n_gene <- round(nrow(Genes)/10)
for(celltype in all_clusters){
  
  DE_sub <- DE %>% 
    filter(cell_type == celltype) %>%
    arrange(padj)
  DE_sub <- DE_sub[1:n_gene,]
  write.table(DE_sub$gene,
              file = paste0(WORK_PATH, "Gene_list_", gsub(" ", "_", celltype), ".txt"),
              sep = "\t", quote = F, col.names = F, row.names = F)
}

# ldcts
LDSC_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Output/"
cluster_file <- gsub(" ", "_", all_clusters)
path_list <- paste0(LDSC_PATH, "Pseudo_Condition_DE/", cluster_file, "_", ",", LDSC_PATH, "Annotation_Gene/ref_")
ldcts <- data.frame(cluster_file, path_list)
write.table(ldcts, row.names = F, col.names = F, quote = F, sep = "\t",
            file = paste0(LDSC_PATH, "Pseudo_Condition_DE/Result/ldcts"))

# DA
n_peak <- round(nrow(Peaks)/10)
for(celltype in all_clusters){
  
  DA_sub <- DA %>% 
    filter(cell_type == celltype) %>%
    arrange(padj)
  DA_sub <- DA_sub[1:n_peak,]
  write.table(DA_sub$peak,
              file = paste0(WORK_PATH, "Peak_list_", gsub(" ", "_", celltype), ".txt"),
              sep = "\t", quote = F, col.names = F, row.names = F)
}

# ldcts
LDSC_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Output/"
cluster_file <- gsub(" ", "_", all_clusters)
path_list <- paste0(LDSC_PATH, "Pseudo_Condition_DA/", cluster_file, "_", ",", LDSC_PATH, "Annotation_Peak2/ref_")
ldcts <- data.frame(cluster_file, path_list)
write.table(ldcts, row.names = F, col.names = F, quote = F, sep = "\t",
            file = paste0(LDSC_PATH, "Pseudo_Condition_DA/Result/ldcts"))


