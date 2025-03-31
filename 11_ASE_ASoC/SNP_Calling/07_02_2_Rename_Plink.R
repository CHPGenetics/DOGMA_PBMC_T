rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)

# Parameters
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/"

for(depth in c(10, 20, 30, 40)){
  
  for(omics in c("RNA", "ATAC", "merged")){
    
    bim <- fread(paste0(WORK_PATH, "Combine/Merged_", omics, "_", depth, ".bim"))
    bim$V2 <- paste0(bim$V1, "_", bim$V4)
    write.table(bim, paste0(WORK_PATH, "Combine/Merged_", omics, "_", depth, ".bim"), 
                col.names = F, row.names = F, sep = "\t", quote = F)
    
  }
  print(depth)
}


