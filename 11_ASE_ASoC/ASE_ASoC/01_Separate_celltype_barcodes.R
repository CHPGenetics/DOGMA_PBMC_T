rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(Seurat)

# Parameters
DATA_PATH <- "/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/01_Barcodes/"

# dogma <- readRDS(paste0(DATA_PATH, "PBMC/output/tcell_annotated_updated.RDS"))
# meta <- dogma@meta.data
# saveRDS(meta, paste0(WORK_PATH, "DOGMA_Meta.rds"))

samples <- fread(paste0(WORK_PATH, "Sample.txt"), header = F)
meta <- readRDS(paste0(WORK_PATH, "DOGMA_Meta.rds"))
meta$celltype_updated <- gsub(") - ", "_", meta$celltype_updated)
meta$celltype_updated <- gsub(" \\(", "_", meta$celltype_updated)
meta$celltype_updated <- gsub(" ", "_", meta$celltype_updated)
meta$celltype_updated <- gsub("\\)", "", meta$celltype_updated)

# write.table(meta$celltype_updated, paste0(WORK_PATH, "CellType.txt"),
#             row.names = F, col.names = F, sep = "\t", quote = F)

for(i in samples$V1){
  
  if(!file.exists(paste0(WORK_PATH, i))){
    
    system(paste0("mkdir ", paste0(WORK_PATH, i)))
  }
  
  for(condi in unique(meta$condition)){
    
    if(!file.exists(paste0(WORK_PATH, i, "/", condi))){
      
      system(paste0("mkdir ", paste0(WORK_PATH, i, "/", condi)))
    }
    for(ctype in unique(meta$celltype_updated)){
      
      meta_sub <- meta %>%
        filter(sample == i) %>%
        filter(condition == condi) %>%
        filter(celltype_updated == ctype)
      if(nrow(meta_sub) > 30){
        
        if(!file.exists(paste0(WORK_PATH, i, "/", condi, "/", ctype))){
          
          system(paste0("mkdir ", paste0(WORK_PATH, i, "/", condi, "/", ctype)))
          system(paste0("mkdir ", paste0(WORK_PATH, i, "/", condi, "/", ctype, "/RNA")))
          system(paste0("mkdir ", paste0(WORK_PATH, i, "/", condi, "/", ctype, "/ATAC")))
        }
        
        write.table(meta_sub$gex_barcode, paste0(WORK_PATH, i, "/", condi, "/", ctype, "/RNA/Barcodes.txt"),
                    row.names = F, col.names = F, sep = "\t", quote = F)
        write.table(meta_sub$atac_barcode, paste0(WORK_PATH, i, "/", condi, "/", ctype, "/ATAC/Barcodes.txt"),
                    row.names = F, col.names = F, sep = "\t", quote = F)
      }
    }
  }
  print(i)
}





