library(data.table)
library(parallel)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(stringr)
library(GenomicRanges)
library(openxlsx)

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

da <- fread('../output/pseudo_bulk_updated/DE_peak_DESeq2_paired.txt')
da <- da[,-c(4,5,9,10,11,13)]
colnames(da)[6:7] <- c('peak', 'nearest_gene')

da_sig <- da

ibd_snps_gr <- readRDS('../output/SNP/ibd_gr.RDS')

merge_gr <- function(peak, gr){
  peak_bp <- peak
  
  peak$chr <- str_split_fixed(peak$peak, '-', 3)[,1]
  peak$start <- str_split_fixed(peak$peak, '-', 3)[,2]
  peak$end <- str_split_fixed(peak$peak, '-', 3)[,3]
  peak_gr <- makeGRangesFromDataFrame(peak, seqnames.field = 'chr', start.field = 'start', end.field = 'end', keep.extra.columns = T)
  
  overlap <- findOverlaps(peak_gr, gr)
  peak_gr_tmp <- peak_gr[queryHits(overlap)]
  gr_tmp <- gr[subjectHits(overlap)]
  
  merge_df <- cbind(as.data.frame(gr_tmp, row.names = 1:length(gr_tmp)), peak = peak_gr_tmp$peak)
  merge_df$snp_pos <- paste(merge_df$seqnames, merge_df$start, sep = '-')
  merge_df <- merge_df[,c(6,12,7,8,9,11)]
  merge_df <- merge_df[!duplicated(merge_df$snp_pos),]
  
  merge_list <- split(merge_df, merge_df$peak)
  
  merge_list_re <- as.data.frame(t(sapply(merge_list, function(df){
    table <- apply(df, 2, function(i){paste0(i, collapse = '|')})
    return(table)
  })))
  merge_list_re$peak <- rownames(merge_list_re)
  
  return(merge_list_re)
}

peak_snp <- merge_gr(da_sig, ibd_snps_gr)
#write.xlsx(peak_snp, '../output/colocalization/peak_snp_updated.xlsx')

da_snp <- merge(da_sig, peak_snp, by = 'peak', all.x = T)

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

p2g <- Reduce(c, p2g_celltype_link)
p2g_df <- as.data.frame(p2g@elementMetadata)
colnames(p2g_df)[1] <- 'Correlation'
p2g_df <- p2g_df[,-c(5,6)]
p2g_df$celltype <- factor(p2g_df$celltype, levels = levels)

p2g_list <- split(p2g_df, p2g_df$celltype)

p2g_list_re <- mclapply(p2g_list, function(df){
  df <- df[,-6]
  df <- df[df$peak %in% da_sig$peak,]
  df_list <- split(df, df$peak)
  
  df_list_re <- as.data.frame(t(sapply(df_list, function(df){
    table <- apply(df, 2, function(i){paste0(i, collapse = '|')})
    return(table)
  })))
  df_list_re$peak <- rownames(df_list_re)
  return(df_list_re)
}, mc.cores = 10)

da_snp$cell_type <- factor(da_snp$cell_type, levels = levels)
da_snp_list <- split(da_snp, da_snp$cell_type)

da_snp_eqta_list <- mclapply(levels[1:20], function(celltype){
  tmp1 <- da_snp_list[[celltype]]
  tmp2 <- p2g_list_re[[celltype]]
  
  tmp <- merge(tmp1, tmp2, all.x = T, by = 'peak')
  return(tmp)
}, mc.cores = 10)
names(da_snp_eqta_list) <- levels

da_snp_eqta <- Reduce(rbind, da_snp_eqta_list)

de <- fread('../output/pseudo_bulk_updated/DE_rna_DESeq2_paired.txt')

de <- de[,-c(4,5,9,10,11)]

de_sig <- de[de$padj < 0.05,]
colnames(de_sig)[2:6] <- paste0('DE_', colnames(de_sig)[2:6])

de_sig$cell_type <- factor(de_sig$cell_type, levels = levels)
de_sig_list <- split(de_sig, de_sig$cell_type)

da_snp_eqta_de_list <- mclapply(levels[1:20], function(celltype){
  tmp1 <- da_snp_eqta_list[[celltype]]
  tmp2 <- de_sig_list[[celltype]]
  tmp2 <- tmp2[,-1]
  tmp2$gene <- tmp2$DE_gene
  tmp2 <- tibble::column_to_rownames(tmp2, var = 'gene')
  
  genes <- str_split(tmp1$gene, '\\|')
  
  eqta_de <- as.data.frame(t(sapply(genes, function(gene){
    shared <- intersect(gene, rownames(tmp2))
    if (length(shared) > 0){
      df <- tmp2[shared,]
      table <- apply(df, 2, function(i){paste0(i, collapse = '|')})
    }
    else{
      table <- rep(NA, 5)
      names(table) <- colnames(tmp2)
    }
    return(table)
  })))
  
  tmp <- cbind(tmp1, eqta_de)
  return(tmp)
}, mc.cores = 10)

da_snp_eqta_de <- Reduce(rbind, da_snp_eqta_de_list)

check_consistency <- function(df){
  peak_log2FoldChange <- df$log2FoldChange
  eqta_gene <- str_split(df$gene, '\\|')
  eqta_correlation <- str_split(df$Correlation, '\\|')
  de_gene <- str_split(df$DE_gene, '\\|')
  de_gene_log2FoldChange <- str_split(df$DE_log2FoldChange, '\\|')
  
  consistent_list <- mclapply(1:nrow(df), function(i){
   if(any(is.na(de_gene[[i]]))){
      gene_consistent <- NA
      gene_length <- 0
    }
    else{
      eqta_correlation_tmp <- as.numeric(eqta_correlation[[i]])[match(de_gene[[i]], eqta_gene[[i]])]
      decision <- peak_log2FoldChange[i]*eqta_correlation_tmp*as.numeric(de_gene_log2FoldChange[[i]]) > 0
      gene_consistent <- de_gene[[i]][decision]
      gene_length <- length(gene_consistent)
      gene_consistent <- paste0(gene_consistent, collapse = '|')
      if (gene_length == 0){
        gene_consistent <- NA
      }
    }
    table <- c(gene_consistent, gene_length)
    return(table)
  }, mc.cores = 10)
  consistent <- as.data.frame(matrix(unlist(consistent_list), ncol = 2, byrow = T))
  colnames(consistent) <- c('gene_consistent', 'gene_consistent_num')
  df_tmp <- cbind(df, consistent)
  return(df_tmp)
}

da_snp_eqta_de_consistent <- check_consistency(da_snp_eqta_de)

write.xlsx(da_snp_eqta_de_consistent, file = '../output/colocalization/not_sig/peak_snp_p2g_gene_sig_consistent_updated.xlsx', overwrite = T)
da_snp_eqta_de_consistent <- da_snp_eqta_de_consistent[da_snp_eqta_de_consistent$gene_consistent_num > 0,]
write.xlsx(da_snp_eqta_de_consistent, file = '../output/colocalization/not_sig/peak_snp_p2g_gene_sig_consistent_only_updated.xlsx', overwrite = T)
da_snp_eqta_de_consistent <- da_snp_eqta_de_consistent[!is.na(da_snp_eqta_de_consistent$snp_pos),]
write.xlsx(da_snp_eqta_de_consistent, file = '../output/colocalization/not_sig/peak_snp_p2g_gene_sig_consistent_only_snp_updated.xlsx', overwrite = T)

check_gene <- function(df, gene){
  genes <- str_split(df$gene_consistent, '\\|')
  tmp <- sapply(genes, function(i){any(i %in% gene)})
  return(tmp)
}

table(unlist(str_split(da_snp_eqta_de_consistent$gene_consistent, '\\|')))[order(table(unlist(str_split(da_snp_eqta_de_consistent$gene_consistent, '\\|'))), decreasing = T)[1:20]]
# PTGER4      TTC33  SATB1-AS1 AC144521.1      SATB1       CREM     GPR183      LNPEP        MYC    RASGRP1 
# 17         14         13         10         10          7          7          7          7          7 
# TNFRSF14      IL2RA  LINC00824     MAP3K8      TAGAP AC009126.1       CAST       ELL2        LIF      RBM17 
# 7          6          6          6          6          5          5          5          5          5 

View(da_snp_eqta_de_consistent[check_gene(da_snp_eqta_de_consistent, 'PTGER4'),])
