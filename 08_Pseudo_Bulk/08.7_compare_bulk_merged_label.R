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

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

ibd_peak <- read.table('../output/SNP/overlap_ibd_peak_snp.txt')

huang <- read.xlsx('../data/SNP/41586_2017_BFnature22969_MOESM2_ESM.xlsx', sheet = 3)
lange <- read.xlsx('../data/SNP/41588_2017_BFng3760_MOESM276_ESM.xlsx')
huang_gene <- huang$Gene
huang_gene <- unlist(str_split(huang_gene, ','))
huang_gene <- huang_gene[!is.na(huang_gene)]
lange_gene <- lange$Implicated.gene
lange_gene <- lange_gene[!is.na(lange_gene)]
lange_gene <- unlist(str_split(lange_gene, ', '))
lange_gene_updated <- c('IL23R', 'ITGA4', 'ITGB8', 'CARD9', 'ITGAL', 'ICAM1', 'TYK2')
lange_gene_updated <- c(lange_gene_updated, lange_gene[-c(1,2,9,14,16,17,22,26,27,28)])
ibd_gene <- union(huang_gene, lange_gene_updated)

rna <- fread('../output/pseudo_bulk/DE_rna_DESeq2_paired_merged.txt')
peak <- fread('../output/pseudo_bulk/DE_peak_DESeq2_paired_merged.txt')

levels = c('CD4+ Naive',
           'CD4+ Regulatory',
           'CD4+ Memory - Th1',
           'CD4+ Memory - Th17',
           'CD4+ Memory - Tfh',
           'CD4+ Memory - Other',
           'CD8+ Naive',
           'CD8+ Regulatory',
           'CD8+ Memory',
           'MAITs', 'Gamma Delta'
)

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
for(i in 1:11){
  celltype <- levels[i]
  print(celltype)
  rna_list[[celltype]] <- result_combine_ps(rna, celltype)
}

bulk <- fread('~/RWorkSpace/Duerr_bulk/RNA/gene/analysis_new/results_Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23.txt')
bulk <- tibble::column_to_rownames(bulk, var = 'V1')
bulk = mutate(bulk, Expression = ifelse(bulk$padj < 0.05, "DEGs", "Not Significant"))
bulk$Expression[is.na(bulk$Expression)] <- "Not Significant"
bulk$Expression[which((bulk$log2FoldChange > 0) & (bulk$padj < 0.05))] <- "Up-regulated"
bulk$Expression[which((bulk$log2FoldChange < 0) & (bulk$padj < 0.05))] <- "Down-regulated"

draw_plot <- function(celltype){
  df1 <- rna_list[[celltype]]
  df2 <- bulk
  df <- merge(df1, df2, all = T, by = 'gene')
  df <- df[order(df$p_val),]
  df$sign <- 'Not Significant'
  df$sign[df$p_val_adj < 0.05] <- 'Significant'
  df$sign[df$p_val_adj < 0.05 & df$padj < 0.05 & df$log2FoldChange*df$avg_logFC > 0] <- 'Replicated'
  df$sign <- factor(df$sign, levels = c('Replicated', 'Significant', 'Not Significant'))
  p <- ggplot(df, aes(y = log2FoldChange, x = avg_logFC))+
    geom_point(aes(color = sign))+
    geom_text_repel(max.overlaps = Inf, aes(label=ifelse(gene %in% ibd_gene & gene %in% df[df$sign == 'Replicated',]$gene,gene,'')))+
    scale_color_manual(values = c("firebrick", alpha("dodgerblue3", 0.7), alpha("lightgrey", 0.1))) +
    theme_bw(base_size = 12) +
    theme(text = element_text(size=16)) +
    theme(legend.position = "bottom") + 
    xlab('Single-cell log2(Fold Change)')+
    ylab('Bulk log2(Fold Change)')+
    labs(col = '', title = celltype, subtitle = paste0(table(df$sign)['Replicated'], ' out of ', table(df$sign)['Replicated']+table(df$sign)['Significant'], ' significant DE genes were replicated in bulk'))+
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

for(i in 1:11){
  celltype <- levels[i]
  print(celltype)
  p <- draw_plot(celltype)
  pdf(paste0('../plots/pseudo_bulk/Bulk_merged_label/RNA/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), draw_plot(celltype))
}

pdf(paste0('../plots/pseudo_bulk/Bulk_merged_label/RNA/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 28, height = 21)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, ncol = 4))
dev.off()

############
peak_list <- list()
for(i in 1:11){
  celltype <- levels[i]
  print(celltype)
  peak_list[[celltype]] <- result_combine_ps(peak, celltype)
}

bulk <- fread('~/RWorkSpace/Duerr_bulk/call_peaks/results_Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23_noHiSeq_annotated.txt')
# bulk <- tibble::column_to_rownames(bulk, var = 'V1')
bulk = mutate(bulk, Expression = ifelse(bulk$padj < 0.05, "DEGs", "Not Significant"))
bulk$Expression[is.na(bulk$Expression)] <- "Not Significant"
bulk$Expression[which((bulk$log2FoldChange > 0) & (bulk$padj < 0.05))] <- "Up-regulated"
bulk$Expression[which((bulk$log2FoldChange < 0) & (bulk$padj < 0.05))] <- "Down-regulated"

draw_plot_peak <- function(celltype){
  df1 <- peak_list[[celltype]]
  df2 <- bulk
  colnames(df2)[10] <- 'gene'
  df <- merge(df1, df2, all = T, by = 'gene')
  df <- df[order(df$p_val),]
  df$sign <- 'Not Significant'
  df$sign[df$p_val_adj < 0.05] <- 'Significant'
  df$sign[df$p_val_adj < 0.05 & df$padj < 0.05 & df$log2FoldChange*df$avg_logFC > 0] <- 'Replicated'
  df$sign <- factor(df$sign, levels = c('Replicated', 'Significant', 'Not Significant'))
  ibd_peak$gene <- ibd_peak$peak_chr_pos
  df <- merge(df, ibd_peak, by = 'gene', all.x = T)
  p <- ggplot(df, aes(y = log2FoldChange, x = avg_logFC))+
    geom_point(aes(color = sign))+
    geom_text_repel(max.overlaps = Inf, aes(label=ifelse(SNP_merged %in% df[df$sign == 'Replicated',]$SNP_merged,SNP_merged,'')))+
    scale_color_manual(values = c("firebrick", alpha("dodgerblue3", 0.7), alpha("lightgrey", 0.1))) +
    theme_bw(base_size = 12) +
    theme(text = element_text(size=16)) +
    theme(legend.position = "bottom") + 
    xlab('Single-cell log2(Fold Change)')+
    ylab('Bulk log2(Fold Change)')+
    labs(col = '', title = celltype, subtitle = paste0(table(df$sign)['Replicated'], ' out of ', table(df$sign)['Replicated']+table(df$sign)['Significant'], ' significant DA peaks were replicated in bulk'))+
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

for(i in 1:11){
  celltype <- levels[i]
  print(celltype)
  p <- draw_plot_peak(celltype)
  pdf(paste0('../plots/pseudo_bulk/Bulk_merged_label/ATAC/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), draw_plot_peak(celltype))
}

pdf(paste0('../plots/pseudo_bulk/Bulk_merged_label/ATAC/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 28, height = 21)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, ncol = 4))
dev.off()

bulk_sig <- bulk[bulk$padj < 0.05,]
table(bulk_sig$name %in% ibd_peak$chr_pos)
# FALSE  TRUE 
# 17002    40 

check_peaks <- function(celltype){
  df <- peak_list[[celltype]]
  df_sig <- df[df$p_val_adj < 0.05,]
  print(table(df_sig$gene %in% ibd_peak$chr_pos))
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

th17 <- rna[rna$cell_type == levels[8],]
