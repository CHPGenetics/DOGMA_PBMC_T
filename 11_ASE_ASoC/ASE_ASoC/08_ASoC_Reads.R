rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(tidyr)

# Parameters
BAM_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/01_Barcodes/"
DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/03_GATK_Allelic_Reads_ATAC/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASoC_Count/"

Samples <- fread("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/Sample_test.txt", header = F)
Samples <- Samples$V1


ASoC_Results <- data.frame()
for(sample in Samples){
  
  for(condition in c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2", 
                     "Act_IL1B_IL23_TGFB", "Act_IL1B_IL23_PGE2_TGFB")){
    
    celltypes <- dir(paste0(BAM_PATH, sample, "/", condition))
    for(celltype in celltypes){
      
      if(file.exists(paste0(BAM_PATH, sample, "/", condition, "/", celltype))){
        
        if(file.exists(paste0(DATA_PATH, condition, "/", sample, "_",
                              condition, "_", celltype, ".gatkAlleleReadCounts.csv"))){
          result <- fread(paste0(DATA_PATH, condition, "/", sample, "_",
                                 condition, "_", celltype, ".gatkAlleleReadCounts.csv")) %>%
            filter(refCount >= 2 & altCount >= 2 & totalCount >= 10) %>%
            # filter(totalCount >= 10) %>%
            mutate(
              p_value = mapply(function(x, n) binom.test(x, n)$p.value, altCount, totalCount),
              p_value_adj = p.adjust(p_value, method = "fdr"),
              sample = sample,
              condition = condition,
              celltype = celltype)
          
          ASoC_Results <- rbind(ASoC_Results, result)
          write.table(ASoC_Results, paste0(WORK_PATH, sample, "_", condition, "_", celltype, "_ASoC_count.csv"), 
                      col.names = T, row.names = F, sep = '\t', quote = F)
        }
        
        print(paste0(condition, ":", sample))
      }
    }
  }
}

write.table(ASoC_Results, paste0(WORK_PATH, "ASoC_counts.txt"), 
            col.names = T, row.names = F, sep = '\t', quote = F)


