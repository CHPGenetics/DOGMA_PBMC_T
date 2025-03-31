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
library(stringr)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

read.dir=function(path){
    files=list.files(path)
    names=str_split_fixed(files, '\\_|\\.', 6)[,5]
    
    file=list()
    for (i in 1:length(files)){
        file[[names[i]]]=fread(paste0(path, files[i]))
    }
    return(file)
}

perturb=read.dir("~/RWorkSpace/GEO_data/DOGMA-seq/asap_reproducibility/CD4_CRISPR_asapseq/output/pseudo_bulk/")

result_combine_ps <- function(de, name){
  data <- de[[name]]
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

############
peak_list <- list()
for(i in 1:5){
  celltype <- names(perturb)[i]
  print(celltype)
  peak_list[[celltype]] <- result_combine_ps(perturb, celltype)
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
    #ibd_peak$gene <- ibd_peak$peak_chr_pos
    #df <- merge(df, ibd_peak, by = 'gene', all.x = T)
    p <- ggplot(df, aes(y = log2FoldChange, x = avg_logFC))+
        geom_point(aes(color = sign))+
        #geom_text_repel(max.overlaps = Inf, aes(label=ifelse(SNP_merged %in% df[df$sign == 'Replicated',]$SNP_merged,SNP_merged,'')))+
        scale_color_manual(values = c(Replicated = "firebrick", Significant = alpha("dodgerblue3", 0.7), `Not Significant` = alpha("lightgrey", 0.1))) +
        theme_bw(base_size = 12) +
        theme(text = element_text(size=16)) +
        theme(legend.position = "bottom") + 
        xlab('Single-cell log2(Fold Change)')+
        ylab('Bulk log2(Fold Change)')+
        labs(col = '', title = celltype, subtitle = paste0(table(df$sign)['Replicated'], ' out of ', table(df$sign)['Replicated']+table(df$sign)['Significant'], ' significant DA peaks were replicated in bulk'))+
        theme(plot.title = element_text(hjust = 0.5))
    return(p)
}

for(i in 1:5){
  celltype <- names(perturb)[i]
  print(celltype)
  p <- draw_plot_peak(celltype)
  pdf(paste0('../plots/pseudo_bulk/Bulk_Perturb/ATAC/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), draw_plot_peak(celltype))
}

pdf(paste0('../plots/pseudo_bulk/Bulk_Perturb/ATAC/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 21, height = 14)
print(wrap_plots(p1, p2, p5, p3, p4, ncol = 3))
dev.off()
