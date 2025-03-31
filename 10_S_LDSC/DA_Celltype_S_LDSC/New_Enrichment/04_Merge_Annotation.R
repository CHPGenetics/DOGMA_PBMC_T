rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(tidyr)

# Parameters
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/03_New_Enrichment/"
SUMM_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Input/"

all_clusters <- readRDS(paste0(SUMM_PATH, "DOGMA_Celltypes.rds"))
Trait <- "IBD"
compare = "01_Celltype_DA/"
cutoff <- "1E-5"

for(i in 1:22){
  
  for(celltype in 1:length(all_clusters)){
    
    if(celltype == 1){
      
      merged_anno <- fread(paste0(WORK_PATH, compare, Trait, "/", cutoff, "/Anno_", 
                                  gsub(" ", "_", all_clusters[celltype]), "_chr", i, ".gz"))
    } else{
      
      anno <- fread(paste0(WORK_PATH, compare, Trait, "/", cutoff, "/Anno_", 
                           gsub(" ", "_", all_clusters[celltype]), "_chr", i, ".gz"))
      merged_anno <- inner_join(merged_anno, anno, by = c("CHR", "SNP", "CM", "BP"))
    }
  }
  
  # # Remove replicate rows
  # merged_anno$sum <- sum(merged_anno[,-1:-4])
  # merged_anno <- merged_anno %>% 
  #   arrange(CHR, BP, desc(sum)) %>%
  #   filter(!duplicated(SNP))
  # merged_anno <- merged_anno[,-"sum"]
  
  # output
  write.table(merged_anno, gzfile(paste0(WORK_PATH, compare, Trait, "/", 
                                         cutoff, "/S_LDSC/Anno_Merged_chr", i, ".annot.gz")), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  print(i)
}
