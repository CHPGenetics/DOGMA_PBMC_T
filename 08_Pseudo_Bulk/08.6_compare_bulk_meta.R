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

rna <- fread('../output/pseudo_bulk_updated/DE_rna_DESeq2_paired.txt')
peak <- fread('../output/pseudo_bulk_updated/DE_peak_DESeq2_paired.txt')

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

rna_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  data <- rna[rna$cell_type == celltype,]
  rna_list[[celltype]] <- result_combine_ps(data)
}

bulk <- fread('~/RWorkSpace/Duerr_bulk/RNA/gene/analysis_new/results_Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23.txt')
bulk <- tibble::column_to_rownames(bulk, var = 'V1')
bulk = mutate(bulk, Expression = ifelse(bulk$padj < 0.05, "DEGs", "Not Significant"))
bulk$Expression[is.na(bulk$Expression)] <- "Not Significant"
bulk$Expression[which((bulk$log2FoldChange > 0) & (bulk$padj < 0.05))] <- "Up-regulated"
bulk$Expression[which((bulk$log2FoldChange < 0) & (bulk$padj < 0.05))] <- "Down-regulated"

meta_analysis <- function(list, celltype, method="DL"){
  lfc_list <- lapply(celltype, function(i){
    return(list[[i]]$log2FoldChange)
  })
  lfc <- Reduce(cbind, lfc_list)
  rownames(lfc) <- list[[1]]$gene
  colnames(lfc) <- celltype
  
  se_list <- lapply(celltype, function(i){
    return(list[[i]]$lfcSE)
  })
  se <- Reduce(cbind, se_list)
  rownames(se) <- list[[1]]$gene
  colnames(se) <- celltype
  
  lfc_no_na <- rownames(lfc)[apply(lfc, 1, function(x){all(!is.na(x))})]
  se_no_na <- rownames(se)[apply(se, 1, function(x){all(!is.na(x))})]
  
  shared <- intersect(lfc_no_na, se_no_na)
  
  lfc_df <- lfc[shared,]
  se_df <- se[shared,]
  
  meta_list <- mclapply(1:nrow(lfc_df), function(i){
    meta_result <- rma(yi = lfc_df[i,], sei = se_df[i,], method=method)
    df <- data.frame(gene = rownames(lfc_df)[i], log2FoldChange = as.numeric(meta_result$beta), lfcSE = meta_result$se, pvalue = meta_result$pval)
    return(df)
  }, mc.cores = 8)
  meta_df <- Reduce(rbind, meta_list)
  
  rownames(meta_df) <- meta_df$gene
  meta_df_df <- meta_df[rownames(lfc),]
  rownames(meta_df_df) <- rownames(lfc)
  meta_df_df$gene <- rownames(meta_df_df)
  meta_df_df$padj <- p.adjust(meta_df_df$pvalue, method = 'BH')
  return(meta_df_df)
}

meta1 <- meta_analysis(rna_list, levels[c(3,4,7,8)], method = 'DL')
meta2 <- meta_analysis(rna_list, levels[c(3,4,7,8)], method = 'FE')
meta3 <- meta_analysis(rna_list, levels[c(4,8)], method = 'DL')
meta4 <- meta_analysis(rna_list, levels[c(4,8)], method = 'FE')

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

p1 <- draw_plot(meta1)
p2 <- draw_plot(meta2)
p3 <- draw_plot(meta3)
p4 <- draw_plot(meta4)

pdf('../plots/pseudo_bulk/meta/RNA/4_random.pdf', width = 7, height = 7)
p1
dev.off()

pdf('../plots/pseudo_bulk/meta/RNA/4_fixed.pdf', width = 7, height = 7)
p2
dev.off()

pdf('../plots/pseudo_bulk/meta/RNA/2_random.pdf', width = 7, height = 7)
p3
dev.off()

pdf('../plots/pseudo_bulk/meta/RNA/2_fixed.pdf', width = 7, height = 7)
p4
dev.off()

pdf('../plots/pseudo_bulk/meta/RNA/meta.pdf', width = 14, height = 14)
wrap_plots(p1,p2,p3,p4, ncol = 2)
dev.off()

############
peak_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  data <- peak[peak$cell_type == celltype,]
  peak_list[[celltype]] <- result_combine_ps(data)
}

bulk <- fread('~/RWorkSpace/Duerr_bulk/call_peaks/results_Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23_noHiSeq_annotated.txt')
# bulk <- tibble::column_to_rownames(bulk, var = 'V1')
bulk = mutate(bulk, Expression = ifelse(bulk$padj < 0.05, "DEGs", "Not Significant"))
bulk$Expression[is.na(bulk$Expression)] <- "Not Significant"
bulk$Expression[which((bulk$log2FoldChange > 0) & (bulk$padj < 0.05))] <- "Up-regulated"
bulk$Expression[which((bulk$log2FoldChange < 0) & (bulk$padj < 0.05))] <- "Down-regulated"
bulk$gene <- bulk$name

meta1 <- meta_analysis(peak_list, levels[c(3,4,7,8)], method = 'DL')
meta2 <- meta_analysis(peak_list, levels[c(3,4,7,8)], method = 'FE')
meta3 <- meta_analysis(peak_list, levels[c(4,8)], method = 'DL')
meta4 <- meta_analysis(peak_list, levels[c(4,8)], method = 'FE')

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
    labs(col = '', subtitle = paste0(table(df$sign)['Replicated'], ' out of ', table(df$sign)['Replicated']+table(df$sign)['Significant'], ' significant DA peaks were replicated in bulk'))+
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

p1 <- draw_plot(meta1)
p2 <- draw_plot(meta2)
p3 <- draw_plot(meta3)
p4 <- draw_plot(meta4)

pdf('../plots/pseudo_bulk/meta/ATAC/4_random.pdf', width = 7, height = 7)
p1
dev.off()

pdf('../plots/pseudo_bulk/meta/ATAC/4_fixed.pdf', width = 7, height = 7)
p2
dev.off()

pdf('../plots/pseudo_bulk/meta/ATAC/2_random.pdf', width = 7, height = 7)
p3
dev.off()

pdf('../plots/pseudo_bulk/meta/ATAC/2_fixed.pdf', width = 7, height = 7)
p4
dev.off()

pdf('../plots/pseudo_bulk/meta/ATAC/meta.pdf', width = 14, height = 14)
wrap_plots(p1,p2,p3,p4, ncol = 2)
dev.off()

