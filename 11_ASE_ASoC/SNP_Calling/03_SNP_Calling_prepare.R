rm(list=ls())
gc()

# Load Packages
library(dplyr)
library(data.table)

# Parameters
DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/02_Combine_BAM/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/"

samples <- fread(paste0(WORK_PATH, "Sample.txt"), header = F)
for(i in 1:nrow(samples)){
  
  sample <- samples[i, 1]
  system(paste0("mkdir ", WORK_PATH, sample))
  
  lst <- c(paste0(sample, "_ATAC,", DATA_PATH, sample, "_ATAC.bam"), 
           paste0(sample, "_RNA,", DATA_PATH, sample, "_GEX.bam"))
  write.table(lst, paste0(WORK_PATH, sample, "/bam.lst"), 
              row.names = F, col.names = F, quote = F, sep = "\t")
}

