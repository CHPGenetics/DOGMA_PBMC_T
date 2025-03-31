rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(dplyr)

# Parameters
pheno = "CD"
DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Data/partition/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Data/LD/"
SUMM_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Input/"
REF_PATH <- "/ix1/wchen/Shiyue/References/LDSC/1000G_Phrase3/1000G_EUR_Phase3_plink/"
PLINK_PATH <- "/ix1/wchen/Shiyue/Biosoft/PLINK/PLINK_1_9/plink"

# SNP ID list
summ <- fread(paste0(SUMM_PATH, pheno, ".txt"))
snp_id <- summ$SNP
write.table(snp_id, paste0(WORK_PATH, pheno, "/snp_id.txt"), row.names = F, col.names = F, quote = F)

# Plink: select SNPs and calculate LD matrix
for(i in 1:22){
  
  system(paste0("mkdir ", WORK_PATH, pheno, "/Chr", i))
  
  KGP <- paste0(REF_PATH, "1000G.EUR.QC.", i)
  bed <- fread(paste0(WORK_PATH, "bed/chr", i, ".txt"), header = F)
  out <- paste0(WORK_PATH, pheno, "/Chr", i, "/chr", i)
  rs_keep <- paste0(WORK_PATH, pheno, "/snp_id.txt")
  select_cmd <- paste0(PLINK_PATH, " --silent --bfile ", 
                       KGP, " --extract ", rs_keep, 
                       " --allow-no-vars --make-bed --recode --out ", out)
  system(select_cmd)
  for(j in 1:nrow(bed)){
    
    bed_partition <- paste0(WORK_PATH, "bed/chr", i, "_", j, ".txt")
    out_ld <- paste0(WORK_PATH, pheno, "/Chr", i, "/chr", i, "_", j)
    ld_cmd <- paste0(PLINK_PATH, " --silent --bfile ",
                     out, " --allow-no-vars --extract range ", bed_partition, 
                     " --make-bed --recode --keep-allele-order --r square --out ", out_ld)
    system(ld_cmd)
  }
  print(i)
}


