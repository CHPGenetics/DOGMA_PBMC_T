library(data.table)
library(parallel)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(venn)
library(ggvenn)
library(ggsci)
library(ggrepel)
library(plyr)
library(doMC)
library(ArchR)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
           'CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
           'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1',
           'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17',
           'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh',
           'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other',
           'CD8+ Naive (Resting)', 'CD8+ Naive (Activated)',
           'CD8+ Regulatory',
           'CD8+ Memory (Resting)', 'CD8+ Memory (Activated)',
           'MAITs (Resting)', 'MAITs (Activated)', 'Gamma Delta'
)

p2g_celltype_link <- list()
for(i in 1:9){
  celltype <- levels[i]
  name <- paste0('0',i,'_',celltype)
  p2g <- readRDS(paste0('../output/ArchR/DOGMA_filtered_multiome/', name, '/p2g_fdr_0.05_loop_replicate_peak_gene.RDS'))
  p2g$Peak2GeneLinks$celltype <- celltype
  p2g_celltype_link[[celltype]] <- p2g$Peak2GeneLinks
}
for(i in 10:20){
  celltype <- levels[i]
  name <- paste0(i,'_',celltype)
  p2g <- readRDS(paste0('../output/ArchR/DOGMA_filtered_multiome/', name, '/p2g_fdr_0.05_loop_replicate_peak_gene.RDS'))
  p2g$Peak2GeneLinks$celltype <- celltype
  p2g_celltype_link[[celltype]] <- p2g$Peak2GeneLinks
}

disease_celltype <- function(disease_snps_gr){
  p2g_tmp <- Reduce(c, p2g_celltype_link)
  overlap <- findOverlaps(disease_snps_gr, p2g_tmp)
  snps_tmp <- disease_snps_gr[queryHits(overlap)]
  p2g_tmp_overlap <- p2g_tmp[subjectHits(overlap)]
  
  snps_df <- cbind(as.data.frame(snps_tmp, row.names = 1:length(snps_tmp)), p2g_tmp_overlap$peak)
  snps_df$snp_pos <- paste(snps_df$seqnames, snps_df$start, sep = '-')
  colnames(snps_df)[9] <- 'peak'
  snps_df <- snps_df[,c(9:10,6:8)]
  snps_df <- snps_df[!duplicated(snps_df$SNP),]
  
  snps_list <- split(snps_df, snps_df$peak)
  
  snps_list_re <- lapply(snps_list, function(df){
    table <- data.frame(peak = names(table(df$peak)),
                        snp_pos = paste(df$snp_pos, collapse = '|'),
                        SNP = paste(df$SNP, collapse = '|'),
                        Disease_merged = paste(df$Disease_merged, collapse = '|'),
                        Group_merged = paste(df$Group_merged, collapse = '|'))
    return(table)
  })
  
  snps_re <- Reduce(rbind, snps_list_re)
  
  p2g_df <- as.data.frame(p2g_tmp@elementMetadata)
  colnames(p2g_df)[1] <- 'Correlation'
  
  overlap_df <- merge(p2g_df, anno, by = 'peak', all.x = T)
  overlap_df_snp <- merge(overlap_df, snps_re, by = 'peak', all.x = T)
  overlap_df_snp <- overlap_df_snp[,c(8,1:7,9:13)]
  overlap_df_snp <- overlap_df_snp[order(overlap_df_snp$FDR),]
  overlap_df_snp$celltype <- factor(overlap_df_snp$celltype, levels = levels)
  overlap_df_snp_df <- Reduce(rbind, split(overlap_df_snp, overlap_df_snp$celltype))
  rownames(overlap_df_snp_df) <- 1:nrow(overlap_df_snp_df)
  
  return(overlap_df_snp_df)
}

disease_celltype_ibd <- function(disease_snps_gr){
  p2g_tmp <- Reduce(c, p2g_celltype_link)
  overlap <- findOverlaps(disease_snps_gr, p2g_tmp)
  snps_tmp <- disease_snps_gr[queryHits(overlap)]
  p2g_tmp_overlap <- p2g_tmp[subjectHits(overlap)]
  
  snps_df <- cbind(as.data.frame(snps_tmp, row.names = 1:length(snps_tmp)), p2g_tmp_overlap$peak)
  snps_df$snp_pos <- paste(snps_df$seqnames, snps_df$start, sep = '-')
  colnames(snps_df)[11] <- 'peak'
  snps_df <- snps_df[,c(11:12,6,8:10)]
  snps_df <- snps_df[!duplicated(snps_df$SNP),]
  
  snps_list <- split(snps_df, snps_df$peak)
  
  snps_list_re <- lapply(snps_list, function(df){
    table <- data.frame(peak = names(table(df$peak)),
                        snp_pos = paste(df$snp_pos, collapse = '|'),
                        SNP = paste(df$SNP, collapse = '|'),
                        trait_uc = paste(df$trait_uc, collapse = '|'),
                        trait_cd = paste(df$trait_cd, collapse = '|'),
                        Disease_merged = paste(df$Disease_merged, collapse = '|'))
    return(table)
  })
  
  snps_re <- Reduce(rbind, snps_list_re)
  
  p2g_df <- as.data.frame(p2g_tmp@elementMetadata)
  colnames(p2g_df)[1] <- 'Correlation'
  
  overlap_df <- merge(p2g_df, anno, by = 'peak', all.x = T)
  overlap_df_snp <- merge(overlap_df, snps_re, by = 'peak', all.x = T)
  overlap_df_snp <- overlap_df_snp[,c(8,1:7,9:14)]
  overlap_df_snp <- overlap_df_snp[order(overlap_df_snp$FDR),]
  overlap_df_snp$celltype <- factor(overlap_df_snp$celltype, levels = levels)
  overlap_df_snp_df <- Reduce(rbind, split(overlap_df_snp, overlap_df_snp$celltype))
  rownames(overlap_df_snp_df) <- 1:nrow(overlap_df_snp_df)
  
  return(overlap_df_snp_df)
}

anno <- fread('../output/pseudo_bulk/peaks.annotation_sorted.txt')[,-1]
anno$peak <- paste(anno$Chr, anno$Start-1, anno$End, sep = '-')
anno <- anno[,c(15,19)]
colnames(anno)[1] <- 'nearest_gene'


disease_snps_gr <- readRDS('../output/SNP/disease_immune_gr.RDS')

table_immune <- disease_celltype(disease_snps_gr)

ibd_snps_gr <- readRDS('../output/SNP/ibd_gr.RDS')

table_ibd <- disease_celltype_ibd(ibd_snps_gr)

library(openxlsx)
write.xlsx(table_immune, file = '../output/colocalization/p2g_snp/p2g_snp_immune.xlsx', colNames = T, rowNames = F)
write.xlsx(table_ibd, file = '../output/colocalization/p2g_snp/p2g_snp_ibd.xlsx', colNames = T, rowNames = F)

table_immune_sub <- table_immune[!is.na(table_immune$snp_pos),]
table_ibd_sub <- table_ibd[!is.na(table_ibd$snp_pos),]
write.xlsx(table_immune_sub, file = '../output/colocalization/p2g_snp/p2g_snp_immune_sub.xlsx', colNames = T, rowNames = F)
write.xlsx(table_ibd_sub, file = '../output/colocalization/p2g_snp/p2g_snp_ibd_sub.xlsx', colNames = T, rowNames = F)
