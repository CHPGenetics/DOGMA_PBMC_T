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

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

ibd_peak <- read.table('../output/SNP/overlap_ibd_peak_snp.txt')

rna <- fread('../output/pseudo_bulk/DE_rna_DESeq2_paired.txt')
peak <- fread('../output/pseudo_bulk/DE_peak_DESeq2_paired.txt')

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

df <- fread('../output/pseudo_bulk/peaks.annotation_sorted.txt')
df$name <- paste(df$Chr, df$Start-1, df$End, sep = '-')
i <- 1
anno <- df

result_combine_ps <- function(de, name){
  data <- de[de$cell_type == name,]
  data <- data[order(data$p_val),]
  results = mutate(data, Expression = ifelse(data$p_val_adj < 0.05, "DEGs", "Not Significant"))
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$avg_logFC > 0) & (results$p_val_adj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$avg_logFC < 0) & (results$p_val_adj < 0.05))] <- "Down-regulated"
  # results$Expression[!results$gene %in% data2[data2$p_val_adj < 0.05,]$gene] <- "Not Significant"
  # p = ggplot(results, aes(avg_log2FC, -log10(p_val_adj))) + geom_point(aes(col = Expression)) +
  #   scale_color_manual(values = c("dodgerblue3", "darkgrey","firebrick")) +
  #   theme_bw(base_size = 12) +
  #   theme(text = element_text(size=16)) +
  #   theme(legend.position = "bottom") + geom_hline(yintercept = -log10(0.05), color = "grey1", lty = 2, lwd = 0.5)+
  #   ylab('-log10(FDR-adjusted P-value)')+
  #   xlab('log2(Fold Change)')+
  #   labs(col = '')+
  #   ggtitle(title)+
  #   theme(plot.title = element_text(hjust = 0.5))
  return(results)
}

rna_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  rna_list[[celltype]] <- result_combine_ps(rna, celltype)
}

############
result_combine_ps_peak <- function(de, name){
  data <- de[de$cell_type == name,]
  print(table(data$gene == anno$name))
  data <- cbind(data, anno[,-c(1:4)])
  data <- data[order(data$p_val),]
  results = mutate(data, Expression = ifelse(data$p_val_adj < 0.05, "DEGs", "Not Significant"))
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$avg_logFC > 0) & (results$p_val_adj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$avg_logFC < 0) & (results$p_val_adj < 0.05))] <- "Down-regulated"
  # results$Expression[!results$gene %in% data2[data2$p_val_adj < 0.05,]$gene] <- "Not Significant"
  # p = ggplot(results, aes(avg_log2FC, -log10(p_val_adj))) + geom_point(aes(col = Expression)) +
  #   scale_color_manual(values = c("dodgerblue3", "darkgrey","firebrick")) +
  #   theme_bw(base_size = 12) +
  #   theme(text = element_text(size=16)) +
  #   theme(legend.position = "bottom") + geom_hline(yintercept = -log10(0.05), color = "grey1", lty = 2, lwd = 0.5)+
  #   ylab('-log10(FDR-adjusted P-value)')+
  #   xlab('log2(Fold Change)')+
  #   labs(col = '')+
  #   ggtitle(title)+
  #   theme(plot.title = element_text(hjust = 0.5))
  return(results)
}

peak_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  peak_list[[celltype]] <- result_combine_ps_peak(peak, celltype)
}

bulk <- fread('~/RWorkSpace/Duerr_bulk/call_peaks/results_Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23_noHiSeq_annotated.txt')
# bulk <- tibble::column_to_rownames(bulk, var = 'V1')
bulk = mutate(bulk, Expression = ifelse(bulk$padj < 0.05, "DEGs", "Not Significant"))
bulk$Expression[is.na(bulk$Expression)] <- "Not Significant"
bulk$Expression[which((bulk$log2FoldChange > 0) & (bulk$padj < 0.05))] <- "Up-regulated"
bulk$Expression[which((bulk$log2FoldChange < 0) & (bulk$padj < 0.05))] <- "Down-regulated"

bulk_sig <- bulk[bulk$padj < 0.05,]
table(bulk_sig$name %in% ibd_peak$peak_chr_pos)
# FALSE  TRUE 
# 17002    40 

check_peaks <- function(celltype){
  df <- peak_list[[celltype]]
  df_sig <- df[df$p_val_adj < 0.05,]
  print(table(df_sig$gene %in% ibd_peak$peak_chr_pos))
}
for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  check_peaks(celltype)
}
# [1] "CD4+ Naive (Resting)"
# 
# FALSE  TRUE 
# 5844    21 
# [1] "CD4+ Naive (Activated)"
# 
# FALSE  TRUE 
# 17615    39 
# [1] "CD4+ Regulatory (Resting)"
# 
# FALSE  TRUE 
# 870     7 
# [1] "CD4+ Regulatory (Activated)"
# 
# FALSE  TRUE 
# 78     1 
# [1] "CD4+ Memory (Resting) - Th1"
# 
# FALSE  TRUE 
# 73754   141 
# [1] "CD4+ Memory (Activated) - Th1"
# 
# FALSE  TRUE 
# 84128   165 
# [1] "CD4+ Memory (Resting) - Th17"
# 
# FALSE  TRUE 
# 56706   114 
# [1] "CD4+ Memory (Activated) - Th17"
# 
# FALSE  TRUE 
# 8003    23 
# [1] "CD4+ Memory (Resting) - Tfh"
# 
# FALSE  TRUE 
# 12560    24 
# [1] "CD4+ Memory (Activated) - Tfh"
# 
# FALSE  TRUE 
# 4715    11 
# [1] "CD4+ Memory (Resting) - Other"
# 
# FALSE  TRUE 
# 35206    71 
# [1] "CD4+ Memory (Activated) - Other"
# 
# FALSE  TRUE 
# 36820    83 
# [1] "CD8+ Naive (Resting)"
# 
# FALSE  TRUE 
# 2609     9 
# [1] "CD8+ Naive (Activated)"
# 
# FALSE  TRUE 
# 9065    26 
# [1] "CD8+ Regulatory"
# 
# FALSE 
# 5 
# [1] "CD8+ Memory (Resting)"
# 
# FALSE  TRUE 
# 34824    75 
# [1] "CD8+ Memory (Activated)"
# 
# FALSE  TRUE 
# 83222   164 
# [1] "MAITs (Resting)"
# 
# FALSE  TRUE 
# 31960    75 
# [1] "MAITs (Activated)"
# 
# FALSE  TRUE 
# 22218    45 
# [1] "Gamma Delta"
# 
# FALSE  TRUE 
# 45669    87 

check_peaks_df <- function(celltype){
  df <- peak_list[[celltype]]
  df_sig <- df[df$p_val_adj < 0.05,]
  ibd_peak$gene <- ibd_peak$peak_chr_pos
  df_sig_snp <- merge(df_sig, ibd_peak[,c(11:15,24)], by = 'gene')
  df_sig_snp <- df_sig_snp[,-c(6:16)]
  return(df_sig_snp)
}
peak_snp_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  peak_snp_list[[celltype]] <- check_peaks_df(celltype)
}
peak_snp <- Reduce(rbind, peak_snp_list)
peak_snp$cell_type <- factor(peak_snp$cell_type, levels = levels)
table(peak_snp$cell_type)
# CD4+ Naive (Resting)          CD4+ Naive (Activated)       CD4+ Regulatory (Resting)     CD4+ Regulatory (Activated) 
# 21                              39                               7                               1 
# CD4+ Memory (Resting) - Th1   CD4+ Memory (Activated) - Th1    CD4+ Memory (Resting) - Th17  CD4+ Memory (Activated) - Th17 
# 141                             165                             114                              23 
# CD4+ Memory (Resting) - Tfh   CD4+ Memory (Activated) - Tfh   CD4+ Memory (Resting) - Other CD4+ Memory (Activated) - Other 
# 24                              11                              71                              83 
# CD8+ Naive (Resting)          CD8+ Naive (Activated)                 CD8+ Regulatory           CD8+ Memory (Resting) 
# 9                              26                               0                              75 
# CD8+ Memory (Activated)                 MAITs (Resting)               MAITs (Activated)                     Gamma Delta 
# 164                              75                              45                              87 

reform_p2g <- function(p2g){
  peak <- metadata(p2g)$peakSet
  gene <- metadata(p2g)$geneSet
  p2g_df <- as.data.frame(p2g)
  
  peak_sub <- peak[p2g_df$idxATAC]
  gene_sub <- gene[p2g_df$idxRNA]
  peak_df <- data.frame(chr = seqnames(peak_sub), start = start(peak_sub), end = end(peak_sub))
  gene_df <- data.frame(chr = seqnames(gene_sub), pos = start(gene_sub), gene = gene_sub$name)
  peak_df$chr_pos <- paste(peak_df$chr, peak_df$start, peak_df$end, sep = '-')
  
  p2g_df$peak_chr <- peak_df$chr
  p2g_df$peak_start <- peak_df$start
  p2g_df$peak_end <- peak_df$end
  p2g_df$peak_chr_pos <- peak_df$chr_pos
  p2g_df$gene_chr <- gene_df$chr
  p2g_df$gene_pos <- gene_df$pos
  p2g_df$gene <- gene_df$gene
  length(unique(p2g_df$peak_chr_pos))#188999
  length(unique(p2g_df$gene))#25736
  p2g_df$distance <- pmin(abs(p2g_df$peak_start - p2g_df$gene_pos), abs(p2g_df$peak_end - p2g_df$gene_pos))
  p2g_df$distance[p2g_df$peak_start <= p2g_df$gene_pos & p2g_df$gene_pos <= p2g_df$peak_end] <- 0
  # p2g_df$gene_merged <- NA
  # p2g_df$gene_pos_merged <- NA
  # p2g_df$gene_num <- NA
  # for (i in 1:nrow(p2g_df)){
  #   tmp <- p2g_df[p2g_df$peak_chr_pos == p2g_df$peak_chr_pos[i],]
  #   p2g_df$gene_merged[i] <- paste(tmp$gene, collapse = '|')
  #   p2g_df$gene_pos_merged[i] <- paste(tmp$gene_pos, collapse = '|')
  #   p2g_df$gene_num[i] <- nrow(tmp)
  # }
  
  # p2g_df_dedup <- p2g_df[!duplicated(p2g_df$peak_chr_pos),]
  # p2g_df_dedup <- p2g_df_dedup[,-c(11:14)]
  
  # merged <- ddply(p2g_df, 'peak_chr_pos', function(tmp){
  #   gene_merged <- paste(tmp$gene, collapse = '|')
  #   distance_merged <- paste(tmp$distance, collapse = '|')
  #   gene_num <- nrow(tmp)
  #   correlation_merged <- paste(tmp$Correlation, collapse = '|')
  #   fdr_merged <- paste(tmp$FDR, collapse = '|')
  #   return(c(gene_merged, distance_merged, gene_num, correlation_merged, fdr_merged))
  # }, .parallel = T)
  # merged <- mclapply(1:nrow(p2g_df_dedup), function(i){
  #   tmp <- p2g_df[p2g_df$peak_chr_pos == p2g_df_dedup$peak_chr_pos[i],]
  #   gene_merged <- paste(tmp$gene, collapse = '|')
  #   gene_pos_merged <- paste(tmp$gene_pos, collapse = '|')
  #   distance_merged <- paste(tmp$distance, collapse = '|')
  #   gene_num <- nrow(tmp)
  #   return(c(gene_merged, gene_pos_merged, distance_merged, gene_num))
  # }, mc.cores = 24)
  # 
  # merged_df <- as.data.frame(matrix(unlist(merged), ncol = 4, byrow = T))
  # colnames(merged_df) <- c('gene_merged', 'gene_pos_merged', 'distance_merged', 'gene_num')
  # colnames(merged)[2:6] <- c('gene_merged', 'distance_merged', 'gene_num', 'correlation_merged', 'fdr_merged')
  p2g_df <- p2g_df[,-c(1:2)]
  return(p2g_df)
}

p2g_celltype_link <- list()
for(i in 1:9){
  celltype <- levels[i]
  name <- paste0('0',i,'_',celltype)
  p2g <- readRDS(paste0('../output/ArchR/DOGMA_filtered_multiome/', name, '/p2g_fdr_0.001.RDS'))
  p2g_celltype_link[[celltype]] <- reform_p2g(p2g)
}
for(i in 10:20){
  celltype <- levels[i]
  name <- paste0(i,'_',celltype)
  p2g <- readRDS(paste0('../output/ArchR/DOGMA_filtered_multiome/', name, '/p2g_fdr_0.001.RDS'))
  p2g_celltype_link[[celltype]] <- reform_p2g(p2g)
}

registerDoMC()
check_peaks_p2g_df <- function(celltype){
  p2g <- p2g_celltype_link[[celltype]]
  peak <- peak_snp_list[[celltype]]
  p2g_sub <- p2g[p2g$peak_chr_pos %in% peak$gene,]
  
  merged <- ddply(p2g_sub, 'peak_chr_pos', function(tmp){
    gene_merged <- paste(tmp$gene, collapse = '|')
    distance_merged <- paste(tmp$distance, collapse = '|')
    gene_num <- nrow(tmp)
    correlation_merged <- paste(tmp$Correlation, collapse = '|')
    fdr_merged <- paste(tmp$FDR, collapse = '|')
    return(c(gene_merged, distance_merged, gene_num, correlation_merged, fdr_merged))
  }, .parallel = T)
  colnames(merged)[2:6] <- c('gene_merged', 'distance_merged', 'gene_num', 'correlation_merged', 'fdr_merged')
  
  merged$gene <- merged$peak_chr_pos
  df <- merge(peak, merged, by = 'gene', all.x = T)
  return(df)
}
peak_snp_p2g_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  peak_snp_p2g_list[[celltype]] <- check_peaks_p2g_df(celltype)
}
peak_snp_p2g_list[[15]] <- peak_snp_p2g_list[[15]][,-c(31:35)]
colnames(peak_snp_p2g_list[[15]]) <- colnames(peak_snp_p2g_list[[1]])
peak_snp_p2g <- Reduce(rbind, peak_snp_p2g_list)
peak_snp_p2g$cell_type <- factor(peak_snp_p2g$cell_type, levels = levels)
table(peak_snp_p2g[!is.na(peak_snp_p2g$gene_num),]$cell_type)
# CD4+ Naive (Resting)          CD4+ Naive (Activated)       CD4+ Regulatory (Resting)     CD4+ Regulatory (Activated) 
# 19                              32                               7                               1 
# CD4+ Memory (Resting) - Th1   CD4+ Memory (Activated) - Th1    CD4+ Memory (Resting) - Th17  CD4+ Memory (Activated) - Th17 
# 83                              91                             103                              23 
# CD4+ Memory (Resting) - Tfh   CD4+ Memory (Activated) - Tfh   CD4+ Memory (Resting) - Other CD4+ Memory (Activated) - Other 
# 24                              11                              64                              67 
# CD8+ Naive (Resting)          CD8+ Naive (Activated)                 CD8+ Regulatory           CD8+ Memory (Resting) 
# 9                              24                               0                              58 
# CD8+ Memory (Activated)                 MAITs (Resting)               MAITs (Activated)                     Gamma Delta 
# 145                              71                              43                              81

p2g_de <- function(genes, table){
  gene <- str_split(genes, '\\|')[[1]]
  table_sub <- tibble::column_to_rownames(table, 'gene')
  table_sub$gene <- rownames(table_sub)
  table_sub <- table_sub[gene,]
  de_p_val_merged <- paste(table_sub$p_val, collapse = '|')
  de_avg_logFC_merged <- paste(table_sub$avg_logFC, collapse = '|')
  # de_pct.1_merged <- paste(table_sub$pct.1, collapse = '|')
  # de_pct.2_merged <- paste(table_sub$pct.2, collapse = '|')
  de_p_val_adj_merged <- paste(table_sub$p_val_adj, collapse = '|')
  tmp <- table_sub$p_val_adj
  if(length(tmp[!is.na(tmp)])>0){
    tmp_min <- min(tmp[!is.na(tmp)])
  }
  else{
    tmp_min <- NA
  }
  de_p_val_adj_min <- tmp_min
  de_Expression_merged <- paste(table_sub$Expression, collapse = '|')
  de_gene_merged <- paste(table_sub$gene, collapse = '|')
  de_gene_sig <- paste(table_sub[table_sub$Expression != 'Not Significant',]$gene, collapse = '|')
  return(c(de_gene_merged, de_avg_logFC_merged, de_p_val_merged, de_p_val_adj_merged, de_p_val_adj_min, de_Expression_merged, de_gene_sig))
}

p2g_de_mclapply <- function(da_p2g, name, table){
  p2g <- da_p2g[[name]]
  merged <- mclapply(1:nrow(p2g), function(i){
    p <- p2g_de(p2g$gene_merged[i], table[[name]])
    return(p)
  }, mc.cores = 4)
  
  merged_df <- as.data.frame(matrix(unlist(merged), ncol = 7, byrow = T))
  colnames(merged_df) <- c('de_gene_merged', 'de_avg_logFC_merged', 'de_p_val_merged', 'de_p_val_adj_merged', 'de_p_val_adj_min', 'de_Expression_merged', 'de_gene_sig')
  rownames(merged_df) <- 1:nrow(merged_df)
  p2g <- cbind(p2g, merged_df)
  return(p2g)
}

peak_snp_p2g_gene_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  peak_snp_p2g_gene_list[[celltype]] <- p2g_de_mclapply(peak_snp_p2g_list, celltype, rna_list)
}
peak_snp_p2g_gene_list[[15]] <- peak_snp_p2g_list[[15]][,c(1:30,1:7)]
colnames(peak_snp_p2g_gene_list[[15]]) <- colnames(peak_snp_p2g_gene_list[[1]])
peak_snp_p2g_gene <- Reduce(rbind, peak_snp_p2g_gene_list)
peak_snp_p2g_gene$cell_type <- factor(peak_snp_p2g_gene$cell_type, levels = levels)

df <- as.matrix(peak_snp_p2g_gene[,c(31:37)])
df[df == 'NA'] <- NA
df[df == ''] <- NA
peak_snp_p2g_gene <- cbind(peak_snp_p2g_gene[,-c(31:37)], as.data.frame(df))

table(peak_snp_p2g_gene[!is.na(peak_snp_p2g_gene$gene_num) & !is.na(peak_snp_p2g_gene$de_gene_sig),]$cell_type)
# CD4+ Naive (Resting)          CD4+ Naive (Activated)       CD4+ Regulatory (Resting)     CD4+ Regulatory (Activated) 
# 12                              24                               2                               0 
# CD4+ Memory (Resting) - Th1   CD4+ Memory (Activated) - Th1    CD4+ Memory (Resting) - Th17  CD4+ Memory (Activated) - Th17 
# 43                              61                              58                              16 
# CD4+ Memory (Resting) - Tfh   CD4+ Memory (Activated) - Tfh   CD4+ Memory (Resting) - Other CD4+ Memory (Activated) - Other 
# 8                               8                              41                              49 
# CD8+ Naive (Resting)          CD8+ Naive (Activated)                 CD8+ Regulatory           CD8+ Memory (Resting) 
# 4                              12                               0                              43 
# CD8+ Memory (Activated)                 MAITs (Resting)               MAITs (Activated)                     Gamma Delta 
# 126                              48                              30                              64 

peak_snp_p2g_gene <- peak_snp_p2g_gene[,-c(7,9:13,15:18,25)]
colnames(peak_snp_p2g_gene)[1] <- 'peak'
df <- peak_snp_p2g_gene[!is.na(peak_snp_p2g_gene$gene_num) & !is.na(peak_snp_p2g_gene$de_gene_sig),]

write.xlsx(peak_snp, file = '../output/colocalization/peak_snp.xlsx', overwrite = T)
write.xlsx(peak_snp_p2g, file = '../output/colocalization/peak_snp_p2g.xlsx', overwrite = T)
write.xlsx(peak_snp_p2g_gene, file = '../output/colocalization/peak_snp_p2g_gene.xlsx', overwrite = T)
write.xlsx(df, file = '../output/colocalization/peak_snp_p2g_gene_sig.xlsx', overwrite = T)

eqta_gene <- str_split(df$gene_merged, '\\|')
eqta_correlation <- str_split(df$correlation_merged, '\\|')
de_gene <- str_split(df$de_gene_merged, '\\|')
de_gene_logFC <- str_split(df$de_avg_logFC_merged, '\\|')
de_gene_adj <- str_split(df$de_p_val_adj_merged, '\\|')

check_consistency <- function(i){
  sign <- df$avg_logFC[i]*as.numeric(eqta_correlation[[i]])*as.numeric(de_gene_logFC[[i]]) > 0
  sig <- as.numeric(de_gene_adj[[i]]) < 0.05
  decision <- sign & sig
  decision[is.na(decision)] <- FALSE
  gene_consistent <- de_gene[[i]][decision]
  gene_length <- length(gene_consistent)
  if (gene_length == 0){
    gene_consistent <- NA
  }
  return(list(gene_consistent, gene_length))
}

gene_consistent <- sapply(1:nrow(df), FUN = function(i){
  paste0(check_consistency(i)[[1]], collapse = '|')
})

gene_consistent_num <- sapply(1:nrow(df), FUN = function(i){
  check_consistency(i)[[2]]
})

df_new <- cbind(df, gene_consistent)
df_new <- cbind(df_new, gene_consistent_num)

write.xlsx(df_new, file = '../output/colocalization/peak_snp_p2g_gene_sig_consistent.xlsx', overwrite = T)
df_new <- df_new[df_new$gene_consistent_num > 0,]
write.xlsx(df_new, file = '../output/colocalization/peak_snp_p2g_gene_sig_consistent_only.xlsx', overwrite = T)
