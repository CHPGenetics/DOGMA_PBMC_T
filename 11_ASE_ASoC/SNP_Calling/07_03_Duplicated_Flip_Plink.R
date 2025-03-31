rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)

# Parameters
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/"

for(depth in c(10, 20, 30, 40)){
  
  for(omics in c("RNA", "ATAC", "merged")){
    
    map <- fread(paste0(WORK_PATH, "Combine/Final/Merged_", omics, "_", depth, ".map"), 
                 col.names = c("chr", "snp_id", "genetic_dist", "pos"), header = F) %>% data.frame()
    map$snp_id <- paste0(map$snp_id, "+", rownames(map))
    
    ped <- fread(paste0(WORK_PATH, "Combine/Final/Merged_", omics, "_", depth, ".ped"), 
                 header = F) %>% data.frame()
    genotype_data <- ped[, 7:ncol(ped)]
    
    for (chr in unique(map$chr)) {
      
      map_chr <- map %>% filter(chr == !!chr)
      # Duplicated SNPs in map file
      dup_pos_chr <- map_chr %>% group_by(pos) %>% filter(n() > 1) %>% pull(pos) %>% unique()
      
      for (dup in dup_pos_chr) {
        # Extract all SNPs in this position
        dup_snps <- map_chr %>% filter(pos == dup)
        
        # Find these SNP index in .ped file
        snp_indices <- map(dup_snps$snp_id, ~which(map$snp_id == .x)) %>% unlist()
        snp_indices_genotype <- c(snp_indices * 2, snp_indices * 2 - 1) %>% sort()
        
        # Extract genotype data
        genotype_at_pos <- genotype_data[, snp_indices_genotype]
        
        # Combine genotype info
        genotype_at_pos_update <- data.frame()
        for (i in 1:nrow(genotype_at_pos)) {
          sample_row <- genotype_at_pos[i, ] %>% t()
          if(all(sample_row == '0')){
            
            genotype_at_pos_update <- rbind(genotype_at_pos_update, c('0', '0'))
            colnames(genotype_at_pos_update) <- colnames(genotype_at_pos)[1:2]
          } else{
            
            non_zero_values <- sample_row[sample_row != 0]
            genotype_at_pos_update <- rbind(genotype_at_pos_update, non_zero_values[1:2])
            colnames(genotype_at_pos_update) <- colnames(genotype_at_pos)[1:2]
          }
        }
        
        genotype_indices_update <- c(snp_indices[1] * 2 - 1, snp_indices[1] * 2)
        genotype_data[, genotype_indices_update] <- genotype_at_pos_update
        genotype_data[,setdiff(snp_indices_genotype, genotype_indices_update)] <- "-"
      }
      print(chr)
    }
    
    map <- map %>% group_by(chr, pos) %>% filter(row_number() == 1)
    map$snp_id <- gsub("\\+.*", "", map$snp_id)
    genotype_data <- genotype_data[, colSums(genotype_data == "-") != nrow(genotype_data)]
    ped <- cbind(ped[, 1:6], genotype_data)
    
    write.table(map, paste0(WORK_PATH, "Combine/Final/Merged_", omics, "_", depth, ".map"), 
                col.names = F, row.names = F, sep = "\t", quote = F)
    write.table(ped, paste0(WORK_PATH, "Combine/Final/Merged_", omics, "_", depth, ".ped"), 
                col.names = F, row.names = F, sep = "\t", quote = F)
    
  }
  print(depth)
}

