library(data.table)
library(parallel)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

rna <- fread('../output/pseudo_bulk/DE_rna_DESeq2_paired_merged.txt')
adt <- fread('../output/pseudo_bulk/DE_adt_DESeq2_paired_merged.txt')
peak <- fread('../output/pseudo_bulk/DE_peak_DESeq2_paired_merged.txt')
library(stringr)
# bed <- str_split_fixed(peak$gene, '-', 3)[1:252371,]
# fwrite(bed, file = '../output/pseudo_bulk/peak.bed', col.names = F, row.names = F, quote = F, sep = '\t')

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

for(i in 1:11){
  celltype <- levels[i]
  print(celltype)
  p <- vol_plot(rna, celltype, 'DE genes')
  pdf(paste0('../plots/pseudo_bulk/RNA_merged/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), vol_plot(rna, celltype, 'DE genes'))
}

pdf(paste0('../plots/pseudo_bulk/RNA_merged/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 28, height = 21)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, ncol = 4))
dev.off()

for(i in 1:11){
  celltype <- levels[i]
  print(celltype)
  p <- vol_plot(adt, celltype, 'DE proteins')
  pdf(paste0('../plots/pseudo_bulk/ADT_merged/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), vol_plot(adt, celltype, 'DE proteins'))
}

pdf(paste0('../plots/pseudo_bulk/ADT_merged/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 28, height = 21)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, ncol = 4))
dev.off()

for(i in 1:11){
  celltype <- levels[i]
  print(celltype)
  p <- vol_plot(peak, celltype, 'DA peaks')
  pdf(paste0('../plots/pseudo_bulk/ATAC_merged/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), vol_plot(peak, celltype, 'DA peaks'))
}

pdf(paste0('../plots/pseudo_bulk/ATAC_merged/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 28, height = 21)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, ncol = 4))
dev.off()

library(openxlsx)
read.marker <- function(path, sheet){
  marker <- read.xlsx(xlsxFile = path, sheet = sheet)
  colnames(marker) <- marker[3,]
  marker <- marker[-c(1:3),]
  marker <- marker %>% mutate(across(contains('mean'), as.numeric)) %>% mutate(across(contains('LFC'), as.numeric)) %>% mutate(across(contains('p_'), as.numeric))
  return(marker)
}

th0.m.marker <- read.marker('../data/41467_2020_15543_MOESM6_ESM.xlsx', sheet = 3)
th0.m.marker.up <- th0.m.marker[th0.m.marker$LFC > 1,]
th0.m.marker.down <- th0.m.marker[th0.m.marker$LFC < -1,]

th0.n.marker <- read.marker('../data/41467_2020_15543_MOESM6_ESM.xlsx', sheet = 1)
th0.n.marker.up <- th0.n.marker[th0.n.marker$LFC > 1,]
th0.n.marker.down <- th0.n.marker[th0.n.marker$LFC < -1,]

th0.m.marker.5d <- read.marker('../data/41467_2020_15543_MOESM6_ESM.xlsx', sheet = 4)
th0.m.marker.5d.up <- th0.m.marker.5d[th0.m.marker.5d$LFC > 1,]
th0.m.marker.5d.down <- th0.m.marker.5d[th0.m.marker.5d$LFC < -1,]

th0.n.marker.5d <- read.marker('../data/41467_2020_15543_MOESM6_ESM.xlsx', sheet = 2)
th0.n.marker.5d.up <- th0.n.marker.5d[th0.n.marker.5d$LFC > 1,]
th0.n.marker.5d.down <- th0.n.marker.5d[th0.n.marker.5d$LFC < -1,]

th0.marker.up <- intersect(intersect(th0.m.marker.up$gene_name, th0.n.marker.up$gene_name), intersect(th0.m.marker.5d.up$gene_name, th0.n.marker.5d.up$gene_name))
th0.marker.down <- intersect(intersect(th0.m.marker.down$gene_name, th0.n.marker.down$gene_name), intersect(th0.m.marker.5d.down$gene_name, th0.n.marker.5d.down$gene_name))

library(ggsci)

vol_plot_label <- function(de, name, title){
  data <- de[de$cell_type == name,]
  results = mutate(data, Expression = ifelse(data$p_val_adj < 0.05, "DEGs", "Not Significant"))
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$avg_logFC > 0) & (results$p_val_adj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$avg_logFC < 0) & (results$p_val_adj < 0.05))] <- "Down-regulated"
  results$Expression[results$Expression != "Not Significant" & results$gene %in% th0.marker.up] <- 'Activation Signature'
  results$Expression[results$Expression != "Not Significant" & results$gene %in% th0.marker.down] <- 'Resting Signature'
  p = ggplot(results, aes(avg_logFC, -log10(p_val_adj))) + geom_point(aes(col = Expression)) +
    scale_color_manual(values = c("Down-regulated" = "#3C5488FF", "Not Significant" = "darkgrey", "Up-regulated" = "#E64B35FF", 'Activation Signature' = "#00A087FF", 'Resting Signature' = "#4DBBD5FF")) +
    theme_bw(base_size = 12) +
    theme(text = element_text(size=16)) +
    theme(legend.position = "bottom") + geom_hline(yintercept = -log10(0.05), color = "grey1", lty = 2, lwd = 0.5)+
    ylab('-log10(FDR-adjusted P-value)')+
    xlab('log2(Fold Change)')+
    labs(col = '', title = name, subtitle = paste0(ifelse(is.na(table(results$Expression)["Up-regulated"]), 0, table(results$Expression)["Up-regulated"])+ifelse(is.na(table(results$Expression)["Down-regulated"]), 0, table(results$Expression)["Down-regulated"]), ' significant ', title, ' were identified\n', ifelse(is.na(table(results$Expression)["Activation Signature"]), 0, table(results$Expression)["Activation Signature"])+ifelse(is.na(table(results$Expression)["Resting Signature"]), 0, table(results$Expression)["Resting Signature"]), ' of them were activation/resting signature genes'))+
    theme(plot.title = element_text(hjust = 0.5))+guides(col = guide_legend(ncol = 3, byrow = T))
  return(p)
}

for(i in 1:11){
  celltype <- levels[i]
  print(celltype)
  p <- vol_plot_label(rna, celltype, 'DE genes')
  pdf(paste0('../plots/pseudo_bulk/RNA_merged_label/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), vol_plot_label(rna, celltype, 'DE genes'))
}

pdf(paste0('../plots/pseudo_bulk/RNA_merged_label/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 28, height = 21)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, ncol = 4))
dev.off()
