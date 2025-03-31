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
library(metafor)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

rna_4 <- fread('../output/pseudo_bulk_updated/DE_rna_DESeq2_paired_4.txt')
peak_4 <- fread('../output/pseudo_bulk_updated/DE_peak_DESeq2_paired_4.txt')
rna_2 <- fread('../output/pseudo_bulk_updated/DE_rna_DESeq2_paired_2.txt')
peak_2 <- fread('../output/pseudo_bulk_updated/DE_peak_DESeq2_paired_2.txt')

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

result_combine_ps <- function(data){
  results = mutate(data, Expression = ifelse(data$padj < 0.05, "DEGs", "Not Significant"))
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$log2FoldChange > 0) & (results$padj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$log2FoldChange < 0) & (results$padj < 0.05))] <- "Down-regulated"
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

bulk <- fread('~/RWorkSpace/Duerr_bulk/RNA/gene/analysis_new/results_Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23.txt')
bulk <- tibble::column_to_rownames(bulk, var = 'V1')
bulk = mutate(bulk, Expression = ifelse(bulk$padj < 0.05, "DEGs", "Not Significant"))
bulk$Expression[is.na(bulk$Expression)] <- "Not Significant"
bulk$Expression[which((bulk$log2FoldChange > 0) & (bulk$padj < 0.05))] <- "Up-regulated"
bulk$Expression[which((bulk$log2FoldChange < 0) & (bulk$padj < 0.05))] <- "Down-regulated"

draw_plot <- function(meta){
  df1 <- result_combine_ps(meta)
  df2 <- bulk
  df <- merge(df1, df2, all = T, by = 'gene')
  df <- df[order(df$pvalue.x),]
  df$sign <- 'Not Significant'
  df$sign[df$padj.x < 0.05] <- 'Significant'
  df$sign[df$padj.x < 0.05 & df$padj.y < 0.05 & df$log2FoldChange.x*df$log2FoldChange.y > 0] <- 'Replicated'
  df$sign <- factor(df$sign, levels = c('Replicated', 'Significant', 'Not Significant'))
  p <- ggplot(df, aes(y = log2FoldChange.y, x = log2FoldChange.x))+
    geom_point(aes(color = sign))+
    scale_color_manual(values = c("firebrick", alpha("dodgerblue3", 0.7), alpha("lightgrey", 0.1))) +
    theme_bw(base_size = 12) +
    theme(text = element_text(size=16)) +
    theme(legend.position = "bottom") + 
    xlab('Single-cell log2(Fold Change)')+
    ylab('Bulk log2(Fold Change)')+
    labs(col = '', subtitle = paste0(table(df$sign)['Replicated'], ' out of ', table(df$sign)['Replicated']+table(df$sign)['Significant'], ' significant DE genes were replicated in bulk'))+
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

p1 <- draw_plot(rna_4)
p2 <- draw_plot(rna_2)

pdf('../plots/pseudo_bulk/merge/RNA/4.pdf', width = 7, height = 7)
p1
dev.off()

pdf('../plots/pseudo_bulk/merge/RNA/2.pdf', width = 7, height = 7)
p2
dev.off()

pdf('../plots/pseudo_bulk/merge/RNA/merge.pdf', width = 14, height = 7)
wrap_plots(p1,p2, ncol = 2)
dev.off()

############
bulk <- fread('~/RWorkSpace/Duerr_bulk/call_peaks/results_Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23_noHiSeq_annotated.txt')
# bulk <- tibble::column_to_rownames(bulk, var = 'V1')
bulk = mutate(bulk, Expression = ifelse(bulk$padj < 0.05, "DEGs", "Not Significant"))
bulk$Expression[is.na(bulk$Expression)] <- "Not Significant"
bulk$Expression[which((bulk$log2FoldChange > 0) & (bulk$padj < 0.05))] <- "Up-regulated"
bulk$Expression[which((bulk$log2FoldChange < 0) & (bulk$padj < 0.05))] <- "Down-regulated"
bulk$gene <- bulk$name

draw_plot <- function(meta){
  df1 <- meta
  df2 <- bulk
  df <- merge(df1, df2, all = T, by = 'gene')
  df <- df[order(df$pvalue.x),]
  df$sign <- 'Not Significant'
  df$sign[df$padj.x < 0.05] <- 'Significant'
  df$sign[df$padj.x < 0.05 & df$padj.y < 0.05 & df$log2FoldChange.x*df$log2FoldChange.y > 0] <- 'Replicated'
  df$sign <- factor(df$sign, levels = c('Replicated', 'Significant', 'Not Significant'))
  p <- ggplot(df, aes(y = log2FoldChange.y, x = log2FoldChange.x))+
    geom_point(aes(color = sign))+
    scale_color_manual(values = c("firebrick", alpha("dodgerblue3", 0.7), alpha("lightgrey", 0.1))) +
    theme_bw(base_size = 12) +
    theme(text = element_text(size=16)) +
    theme(legend.position = "bottom") + 
    xlab('Single-cell log2(Fold Change)')+
    ylab('Bulk log2(Fold Change)')+
    labs(col = '', subtitle = paste0(table(df$sign)['Replicated'], ' out of ', table(df$sign)['Replicated']+table(df$sign)['Significant'], ' significant DA peaks were replicated in bulk'))+
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

p1 <- draw_plot(peak_4)
p2 <- draw_plot(peak_2)

pdf('../plots/pseudo_bulk/merge/ATAC/4.pdf', width = 7, height = 7)
p1
dev.off()

pdf('../plots/pseudo_bulk/merge/ATAC/2.pdf', width = 7, height = 7)
p2
dev.off()

pdf('../plots/pseudo_bulk/merge/ATAC/merge.pdf', width = 14, height = 7)
wrap_plots(p1,p2, ncol = 2)
dev.off()
