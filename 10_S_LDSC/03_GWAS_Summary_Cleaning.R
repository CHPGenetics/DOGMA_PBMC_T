rm(list=ls())
gc()

# Load Packages
library(data.table)
library(plyr)
library(dplyr)
library(rtracklayer)
library(liftOver)
library(GenomicRanges)

# Parameters
DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Data/GWAS_Summary_Raw/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Input/"
REF_PATH <- "/ix1/wchen/Shiyue/References/LDSC/"

# Munge Function
# load reference
hm3_snplist <- fread(paste0(REF_PATH, "w_hm3.snplist"))
hm3_snp_info <- llply(c(1: 22), function(chr){
  
  ldsc_chr <- fread(paste0(REF_PATH, "eur_w_ld_chr/", chr, ".l2.ldscore.gz"))
}) %>% do.call("rbind", .) %>% select((CHR: BP))
hm3_snp_info_join <- left_join(hm3_snplist, hm3_snp_info, by = "SNP")

# munge function
munge.hm3 <- function(data_ldsc){
  
  data_ldsc_ref <- inner_join(hm3_snplist, data_ldsc, by = c("SNP"))
  ## delete allele
  na_row1 <- which(data_ldsc_ref$A1.x != data_ldsc_ref$A2.y & data_ldsc_ref$A1.x != data_ldsc_ref$A1.y)
  na_row2 <- which(data_ldsc_ref$A2.x != data_ldsc_ref$A2.y & data_ldsc_ref$A2.x != data_ldsc_ref$A1.y)
  na_row <- c(na_row1, na_row2) %>% unique()
  if(length(na_row) != 0){
    
    data_ldsc_ref$N[na_row] <- NA
    data_ldsc_ref$Z[na_row] <- NA
  }
  ## flip allele
  fl_row <- which(data_ldsc_ref$A1.x == data_ldsc_ref$A2.y & data_ldsc_ref$A2.x == data_ldsc_ref$A1.y)
  if(length(fl_row) != 0){
    
    data_ldsc_ref$Z[fl_row] <- -1 * data_ldsc_ref$Z[fl_row]
  }
  data_ldsc_ref_f <- data.frame(SNP = data_ldsc_ref$SNP, 
                                CHR = data_ldsc_ref$CHR,
                                BP = data_ldsc_ref$BP,
                                A1 = data_ldsc_ref$A1.x,
                                A2 = data_ldsc_ref$A2.x,
                                Z = data_ldsc_ref$Z,
                                P = data_ldsc_ref$P,
                                N = data_ldsc_ref$N)
  
  return(data_ldsc_ref_f)
}

# Process Summary
## 01_IBD
IBD <- fread(paste0(DATA_PATH, "NG_IBD_3760/ibd_build37_59957_20161107.txt"))
IBD <- IBD %>%
  separate(col = MarkerName, into = c("CHR", "BP"), sep = ":") %>%
  mutate(BP = gsub("_.*", "", BP))
IBD <- data.frame(CHR = as.numeric(IBD$CHR),
                  BP = as.numeric(IBD$BP),
                  N = 59957,
                  Z = as.numeric(IBD$Effect/IBD$StdErr),
                  P = as.numeric(IBD$P.value),
                  A1 = toupper(IBD$Allele1),
                  A2 = toupper(IBD$Allele2))
IBD <- inner_join(hm3_snp_info, IBD, by = c("CHR", "BP"))
IBD <- IBD[!duplicated(IBD$SNP),]
IBD <- munge.hm3(IBD)
write.table(IBD, row.names = F, quote = F, na = "NA", sep = "\t",
            file = paste0(WORK_PATH, "IBD.txt"))

## 02_CD
CD <- fread(paste0(DATA_PATH, "NG_IBD_3760/cd_build37_40266_20161107.txt"))
CD <- CD %>%
  separate(col = MarkerName, into = c("CHR", "BP"), sep = ":") %>%
  mutate(BP = gsub("_.*", "", BP))
CD <- data.frame(CHR = as.numeric(CD$CHR),
                  BP = as.numeric(CD$BP),
                  N = 40266,
                  Z = as.numeric(CD$Effect/CD$StdErr),
                  P = as.numeric(CD$P.value),
                  A1 = toupper(CD$Allele1),
                  A2 = toupper(CD$Allele2))
CD <- inner_join(hm3_snp_info, CD, by = c("CHR", "BP"))
CD <- CD[!duplicated(CD$SNP),]
CD <- munge.hm3(CD)
write.table(CD, row.names = F, quote = F, na = "NA", sep = "\t",
            file = paste0(WORK_PATH, "CD.txt"))

## 03_UC
UC <- fread(paste0(DATA_PATH, "NG_IBD_3760/uc_build37_45975_20161107.txt"))
UC <- UC %>%
  separate(col = MarkerName, into = c("CHR", "BP"), sep = ":") %>%
  mutate(BP = gsub("_.*", "", BP))
UC <- data.frame(CHR = as.numeric(UC$CHR),
                  BP = as.numeric(UC$BP),
                  N = 45975,
                  Z = as.numeric(UC$Effect/UC$StdErr),
                  P = as.numeric(UC$P.value),
                  A1 = toupper(UC$Allele1),
                  A2 = toupper(UC$Allele2))
UC <- inner_join(hm3_snp_info, UC, by = c("CHR", "BP"))
UC <- UC[!duplicated(UC$SNP),]
UC <- munge.hm3(UC)
write.table(UC, row.names = F, quote = F, na = "NA", sep = "\t",
            file = paste0(WORK_PATH, "UC.txt"))





