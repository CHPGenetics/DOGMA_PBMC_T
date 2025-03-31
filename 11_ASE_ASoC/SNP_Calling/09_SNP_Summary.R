rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)

# Parameters
DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/Combine/Final/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/04_SNP_Summary/"

# Overall
summary <- data.frame()
for(depth in c(10, 20, 30, 40)){
  
  for(omics in c("RNA", "ATAC", "merged")){
    
    for(kgp in c("", "_KGP", "_Geno", "_Geno_KGP")){
      
      if(kgp == ""){
        
        maf <- ""
        bim <- fread(paste0(DATA_PATH, omics, "_", depth, ".bim"))
        n <- nrow(bim)
        info <- c(Resource = omics, Depth = depth, KGP = kgp, MAF = maf, N = n)
        
      } 
      if(kgp %in% c("_Geno", "_Geno_KGP")){
        
        maf <- ""
        bim <- fread(paste0(DATA_PATH, omics, "_", depth, kgp, ".bim"))
        n <- nrow(bim)
        info <- c(Resource = omics, Depth = depth, KGP = kgp, MAF = maf, N = n)
        
      } 
      if(kgp == "_KGP"){
        
        info <- data.frame()
        for(maf in c("", "_0.01", "_0.05")){
          
          bim <- fread(paste0(DATA_PATH, omics, "_", depth, kgp, maf, ".bim"))
          n <- nrow(bim)
          info_s <- c(Resource = omics, Depth = depth, KGP = kgp, MAF = maf, N = n)
          info <- rbind(info, info_s)
          colnames(info) <- c("Resource", "Depth", "KGP", "MAF", "N")
        }
      }
      
      summary <- rbind(summary, info)
      colnames(summary) <- c("Resource", "Depth", "KGP", "MAF", "N")
    }
  }
  print(depth)
}

summary <- summary %>%
  mutate(KGP = ifelse(KGP == "_KGP", "KGP", 
                      ifelse(KGP == "_Geno", "Geno 0.5", 
                             ifelse(KGP == "_Geno_KGP", "Geno KGP", "All")))) %>%
  mutate(MAF = ifelse(MAF == "_0.01", 0.01, 
                      ifelse(MAF == "_0.05", 0.05, 0)))

write.table(summary, paste0(WORK_PATH, "SNP_Number_overall.txt"), 
            col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(summary, paste0(WORK_PATH, "SNP_Number_overall.csv"))


# Calculate the proportion
proportion <- data.frame()
summary_sub <- summary %>%
  filter(KGP %in% c("KGP", "All"))
for(depth in c(10, 20, 30, 40)){
  
  for(omics in c("RNA", "ATAC", "merged")){
    
    kgp <- summary_sub$N[which(summary_sub$Resource == omics &
                             summary_sub$Depth == depth & summary_sub$KGP == "KGP" & 
                               summary_sub$MAF == 0)] %>%
      as.numeric()
    overall <- summary_sub$N[which(summary_sub$Resource == omics &
                                 summary_sub$Depth == depth & summary_sub$KGP == "All")] %>%
      as.numeric()
    prop <- kgp/overall
    info <- c(Resource = omics, Depth = depth, Type = "All", KGP_n = kgp, Prop = prop, 
              Proportion = paste0(round(prop * 100, 2), "%"))
    proportion <- rbind(proportion, info)
    colnames(proportion) <- c("Resource", "Depth", "Type", "KGP_N", "Prop", "Proportion")
  }
}

summary_sub <- summary %>%
  filter(KGP %in% c("Geno KGP", "Geno 0.5"))
for(depth in c(10, 20, 30, 40)){
  
  for(omics in c("RNA", "ATAC", "merged")){
    
    kgp <- summary_sub$N[which(summary_sub$Resource == omics &
                                 summary_sub$Depth == depth & summary_sub$KGP == "Geno KGP" & 
                                 summary_sub$MAF == 0)] %>%
      as.numeric()
    overall <- summary_sub$N[which(summary_sub$Resource == omics &
                                     summary_sub$Depth == depth & summary_sub$KGP == "Geno 0.5")] %>%
      as.numeric()
    prop <- kgp/overall
    info <- c(Resource = omics, Depth = depth, Type = "Geno", KGP_n = kgp, Prop = prop, 
              Proportion = paste0(round(prop * 100, 2), "%"))
    proportion <- rbind(proportion, info)
    colnames(proportion) <- c("Resource", "Depth", "Type", "KGP_N", "Prop", "Proportion")
  }
}

write.table(proportion, paste0(WORK_PATH, "SNP_Proportion_overall.txt"), 
            col.names = T, row.names = F, sep = "\t", quote = F)

proportion <- proportion %>% filter(Type == "Geno")
write.csv(proportion, paste0(WORK_PATH, "SNP_Proportion_overall.csv"), row.names = F)
