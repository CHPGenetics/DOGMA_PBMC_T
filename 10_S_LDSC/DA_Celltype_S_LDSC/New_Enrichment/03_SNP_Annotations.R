rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(tidyr)
library(readxl)

# Parameters
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/03_New_Enrichment/"
SUMM_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Input/"
REF_PATH <- "/ix1/wchen/Shiyue/References/LDSC/1000G_Phrase3/1000G_EUR_Phase3_plink/"
DATA_PATH <- "/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/pseudo_bulk_updated/"
PEAK_REF_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Input/"

pheno <- "IBD"

# DA P value transformation
Peaks <- read.table(paste0(PEAK_REF_PATH, "Peak3_coord_all.txt"), header = T)
all_clusters <- readRDS(paste0(SUMM_PATH, "DOGMA_Celltypes.rds"))

colnames(Peaks) <- c("Peak", "CHR", "Start", "End")
for(celltype in all_clusters){
  
  DA_sub <- fread(paste0(PEAK_REF_PATH, "Pseudo_Celltype_DA/", "DA_", gsub(" ", "_", celltype), ".txt")) %>% 
    mutate(cell_type = celltype) %>%
    filter(Peak %in% Peaks$Peak) %>%
    arrange(p_val)
  
  DA_sub$p_val[which(DA_sub$p_val == 0) ] <- min(DA_sub$p_val[DA_sub$p_val != 0])
  
  # Min/Max Trans
  DA_sub <- DA_sub %>% filter(!is.na(p_val)) %>% 
    mutate(p_val_trans = -2*log(p_val)) %>%
    mutate(p_val_trans = (p_val_trans - min(p_val_trans))/(max(p_val_trans) - min(p_val_trans)))
  
  # DA_sub <- separate(DA_sub, peak, into = c("CHR", "Start", "End"), sep = "-", remove = FALSE)
  # DA_sub$CHR <- gsub("chr", "", DA_sub$CHR) %>% as.numeric()
  DA_sub <- left_join(DA_sub, Peaks, by = "Peak")
  
  DA_sub <- DA_sub %>% 
    mutate(Start = as.numeric(Start)) %>%
    mutate(End = as.numeric(End)) %>%
    select(Peak, CHR, Start, End, p_val_trans)
  write.table(DA_sub,
              file = paste0(WORK_PATH, "00_", pheno, "_SuSIE_Celltype/Peak_list_", gsub(" ", "_", celltype), ".txt"),
              sep = "\t", quote = F, col.names = T, row.names = F)
  print(celltype)
}

# Create annotate.gz
# Add fine-mapping PIP
compare = "01_Celltype_DA/"
cutoff <- "1E-5"
windowsize <- 100000
PIP <- read.table(paste0(WORK_PATH, "00_", pheno, "_SuSIE/SuSIE_PIP_", cutoff, ".txt"), header = T)
PIP <- PIP %>% arrange(CHR, BP, desc(PIP)) %>%
  filter(!duplicated(SNP, CHR, BP))

for(celltype in all_clusters){
  
  DA_sub <- fread(paste0(WORK_PATH, "00_", pheno, "_SuSIE_Celltype/Peak_list_", gsub(" ", "_", celltype), ".txt"))
  for(i in 1:22){
    
    bim <- fread(paste0(REF_PATH, "1000G.EUR.QC.", i, ".bim"))
    bim <- select(bim, 1:4)
    
    # Calculate the weight
    DA_sub_i <- DA_sub %>% filter(CHR == i)
    PIP_i <- PIP %>% filter(CHR == i)
    
    if(nrow(PIP_i) > 0){
      
      PIP_i <- PIP_i %>%
        mutate(P_trans = -2*log(P),
               PIP_trans = 2*log(PIP + 1)) %>%
        mutate(P_trans = (P_trans - min(P_trans))/(max(P_trans) - min(P_trans)),
               PIP_trans = (PIP_trans - min(PIP_trans))/(max(PIP_trans) - min(PIP_trans))) %>%
        mutate(Weight = sqrt((P_trans + PIP_trans)/2))
      
      SNP_weights <- PIP_i %>%
        # For each SNP
        rowwise() %>%
        do({
          SNP_row <- .
          matching_peaks <- DA_sub_i %>%
            filter(SNP_row$BP >= (Start - windowsize) & SNP_row$BP <= (End + windowsize)) %>%
            arrange(desc(p_val_trans)) %>%
            slice(1) %>%
            mutate(weight = sqrt(p_val_trans * SNP_row$Weight), SNP = SNP_row$SNP) %>% # add SNP column
            select(SNP, Peak, weight) # 确保包含SNP列
          
          if (nrow(matching_peaks) == 0) {
            data.frame(SNP = SNP_row$SNP, Peak = NA, weight = 0)
          } else {
            matching_peaks
          }
        }) %>%
        unnest()
      
      # Combine
      SNP_weights <- SNP_weights[,-2]
      
      bim_anno <- left_join(bim, SNP_weights, by = c("V2" = "SNP")) %>%
        mutate(weight = replace_na(weight, 0))
      colnames(bim_anno) <- c("CHR", "SNP", "CM", "BP", gsub(" ", "_", celltype))
      
      # Output
      write.table(bim_anno, gzfile(paste0(WORK_PATH, compare, pheno, "/", cutoff, "/Anno_", gsub(" ", "_", celltype), 
                                          "_chr", i, ".gz")), sep = "\t", row.names = FALSE, 
                  col.names = TRUE, quote = FALSE)
    } else {
      
      SUMM_i <- fread(paste0(SUMM_PATH, pheno, ".txt")) %>%
        filter(CHR == 1)
      SNP_weights <- data.frame(SNP = SUMM_i$SNP,
                                weight = 0)
      bim_anno <- left_join(bim, SNP_weights, by = c("V2" = "SNP")) %>%
        mutate(weight = replace_na(weight, 0))
      colnames(bim_anno) <- c("CHR", "SNP", "CM", "BP", gsub(" ", "_", celltype))
      
      write.table(bim_anno, gzfile(paste0(WORK_PATH, compare, pheno, "/", cutoff, "/Anno_", gsub(" ", "_", celltype), 
                                          "_chr", i, ".gz")), sep = "\t", row.names = FALSE, 
                  col.names = TRUE, quote = FALSE)
    }
  }
}
