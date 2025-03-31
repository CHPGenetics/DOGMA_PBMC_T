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


check_overlap_5 <- function(table1, table2){
  df1 <- table2$`CD4+ Memory - Th1`
  df2 <- table2$`CD4+ Memory - Th1`
  df3 <- table2$`CD4+ Memory - Th17`
  df4 <- table2$`CD4+ Memory - Th17`
  df5 <- table2$`CD4+ Regulatory`
  x <- list(
    `Bulk RNA-seq` = table1[table1$Expression!='Not Significant',]$gene, 
    `CD4+ Memory - Th1` = df1[df1$Expression!='Not Significant',]$gene,
    `CD4+ Memory - Th1` = df2[df2$Expression!='Not Significant',]$gene,
    `CD4+ Memory - Th17` = df3[df3$Expression!='Not Significant',]$gene,
    `CD4+ Memory - Th17` = df4[df4$Expression!='Not Significant',]$gene,
    `CD4+ Regulatory` = df5[df5$Expression!='Not Significant',]$gene
  )
  # p1 <- ggvenn(
  #   x, 
  #   fill_color = pal_npg()(6),
  #   stroke_size = 0.5, set_name_size = 4
  # )
  p1 <- venn(x, ilab=TRUE, zcolor = "style", ggplot = T)
  return(p1)
}

p <- check_overlap_5(bulk, rna_list)

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
    geom_text_repel(max.overlaps = Inf, aes(label=ifelse(gene %in% df[df$sign == 'Replicated',]$gene[1:10],gene,'')))+
    scale_color_manual(values = c("firebrick", alpha("dodgerblue3", 0.7), alpha("darkgrey", 0.3))) +
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
  pdf(paste0('../plots/pseudo_bulk/Bulk_merged/RNA/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), draw_plot(celltype))
}

pdf(paste0('../plots/pseudo_bulk/Bulk_merged/RNA/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 28, height = 21)
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
  p <- ggplot(df, aes(y = log2FoldChange, x = avg_logFC))+
    geom_point(aes(color = sign))+
    scale_color_manual(values = c("firebrick", alpha("dodgerblue3", 0.7), alpha("darkgrey", 0.3))) +
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
  pdf(paste0('../plots/pseudo_bulk/Bulk_merged/ATAC/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), draw_plot_peak(celltype))
}

pdf(paste0('../plots/pseudo_bulk/Bulk_merged/ATAC/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 28, height = 21)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, ncol = 4))
dev.off()
