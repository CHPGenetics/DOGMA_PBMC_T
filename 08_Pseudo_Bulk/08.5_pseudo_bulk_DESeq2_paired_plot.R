library(data.table)
library(parallel)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

rna <- fread('../output/pseudo_bulk/DE_rna_DESeq2_paired.txt')
adt <- fread('../output/pseudo_bulk/DE_adt_DESeq2_paired.txt')
peak <- fread('../output/pseudo_bulk/DE_peak_DESeq2_paired.txt')
library(stringr)
bed <- str_split_fixed(peak$gene, '-', 3)[1:252371,]
fwrite(bed, file = '../output/pseudo_bulk/peak.bed', col.names = F, row.names = F, quote = F, sep = '\t')

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

vol_plot <- function(de, name, title){
  data <- de[de$cell_type == name,]
  results = mutate(data, Expression = ifelse(data$p_val_adj < 0.05, "DEGs", "Not Significant"))
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$avg_logFC > 0) & (results$p_val_adj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$avg_logFC < 0) & (results$p_val_adj < 0.05))] <- "Down-regulated"
  p = ggplot(results, aes(avg_logFC, -log10(p_val_adj))) + geom_point(aes(col = Expression)) +
    scale_color_manual(values = c("Down-regulated" = "dodgerblue3", "Not Significant" = "darkgrey", "Up-regulated" = "firebrick")) +
    theme_bw(base_size = 12) +
    theme(text = element_text(size=16)) +
    theme(legend.position = "bottom") + geom_hline(yintercept = -log10(0.05), color = "grey1", lty = 2, lwd = 0.5)+
    ylab('-log10(FDR-adjusted P-value)')+
    xlab('log2(Fold Change)')+
    labs(col = '', title = name, subtitle = paste0(ifelse(is.na(table(results$Expression)["Up-regulated"]), 0, table(results$Expression)["Up-regulated"])+ifelse(is.na(table(results$Expression)["Down-regulated"]), 0, table(results$Expression)["Down-regulated"]), ' significant ', title, ' were identified'))+
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  p <- vol_plot(rna, celltype, 'DE genes')
  pdf(paste0('../plots/pseudo_bulk/RNA/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), vol_plot(rna, celltype, 'DE genes'))
}

pdf(paste0('../plots/pseudo_bulk/RNA/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 35, height = 28)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, ncol = 5))
dev.off()

for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  p <- vol_plot(adt, celltype, 'DE proteins')
  pdf(paste0('../plots/pseudo_bulk/ADT/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), vol_plot(adt, celltype, 'DE proteins'))
}

pdf(paste0('../plots/pseudo_bulk/ADT/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 35, height = 28)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, ncol = 5))
dev.off()

for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  p <- vol_plot(peak, celltype, 'DA peaks')
  pdf(paste0('../plots/pseudo_bulk/ATAC/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), vol_plot(peak, celltype, 'DA peaks'))
}

pdf(paste0('../plots/pseudo_bulk/ATAC/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 35, height = 28)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, ncol = 5))
dev.off()

th17 <- peak[peak$cell_type == 'CD4+ Memory (Activated) - Th17',]
df <- fread('../output/pseudo_bulk/peaks.annotation_sorted.txt')
df$name <- paste(df$Chr, df$Start-1, df$End, sep = '-')
i <- 1
anno <- df
while (i < 20){
  anno <- rbind(anno, df)
  i <- i + 1
}
table(peak$gene == anno$name)
peak_anno <- cbind(peak, anno[,-c(1:4)])
fwrite(peak_anno, file = '../output/pseudo_bulk/DE_peak_DESeq2_paired_anno.txt', col.names = T, row.names = F, quote = F)
