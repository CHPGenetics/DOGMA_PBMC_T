library(data.table)
library(parallel)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(ggrepel)
library(openxlsx)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

rna <- fread('../output/pseudo_bulk_updated/DE_rna_DESeq2_paired.txt')
adt <- fread('../output/pseudo_bulk_updated/DE_adt_DESeq2_paired.txt')
peak <- fread('../output/pseudo_bulk_updated/DE_peak_DESeq2_paired.txt')
library(stringr)

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
  results = mutate(data, Expression = ifelse(data$padj < 0.05, "DEGs", "Not Significant"))
  results = results[order(results$pvalue),]
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$log2FoldChange > 0) & (results$padj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$log2FoldChange < 0) & (results$padj < 0.05))] <- "Down-regulated"
  p = ggplot(results, aes(log2FoldChange, -log10(padj))) + geom_point(aes(col = Expression)) +
    geom_text_repel(max.overlaps = Inf, min.segment.length = -1, aes(label=ifelse(gene %in% results[results$Expression != "Not Significant",]$gene[1:20],gene,'')))+
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

ipa <- function(de, name){
  data <- de[de$cell_type == name,]
  results = mutate(data, Expression = ifelse(data$padj < 0.05, "DEGs", "Not Significant"))
  results = results[order(results$pvalue),]
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$log2FoldChange > 0) & (results$padj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$log2FoldChange < 0) & (results$padj < 0.05))] <- "Down-regulated"
  return(results)
}

for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  p <- vol_plot(rna, celltype, 'DE genes')
  pdf(paste0('../plots/pseudo_bulk/RNA_updated/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), vol_plot(rna, celltype, 'DE genes'))
}

pdf(paste0('../plots/pseudo_bulk/RNA_updated/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 35, height = 28)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, ncol = 5))
dev.off()

for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  results <- ipa(rna, celltype)
  results = results[results$Expression != "Not Significant",]
  write.xlsx(results, file = paste0('../plots/pseudo_bulk/RNA_updated/IPA/', celltype, '.xlsx'), colNames = T, rowNames = F)
}

for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  p <- vol_plot(adt, celltype, 'DE proteins')
  pdf(paste0('../plots/pseudo_bulk/ADT_updated/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), vol_plot(adt, celltype, 'DE proteins'))
}

pdf(paste0('../plots/pseudo_bulk/ADT_updated/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 35, height = 28)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, ncol = 5))
dev.off()

for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  p <- vol_plot(peak, celltype, 'DA peaks')
  pdf(paste0('../plots/pseudo_bulk/ATAC_updated/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), vol_plot(peak, celltype, 'DA peaks'))
}

pdf(paste0('../plots/pseudo_bulk/ATAC_updated/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 35, height = 28)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, ncol = 5))
dev.off()

peak_p2g = read.xlsx('../output/colocalization/not_sig/peak_snp_p2g_gene_sig_consistent_updated.xlsx')

vol_plot_new <- function(de, name, title){
  data <- de[de$cell_type == name,]
  results = mutate(data, Expression = ifelse(data$padj < 0.05, "DEGs", "Not Significant"))
  results = results[order(results$pvalue),]
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$log2FoldChange > 0) & (results$padj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$log2FoldChange < 0) & (results$padj < 0.05))] <- "Down-regulated"
  
  tmp1 = str_split(results$FDR, '\\|')
  index = lapply(tmp1, function(i){which.min(as.numeric(i))})
  tmp2 = str_split(results$gene, '\\|')

  top_linked_gene_list = lapply(1:length(index), function(i){
    tmp3 = index[[i]]
    if(length(tmp3) == 0){
      tmp4 = NA
    }
    else{
      tmp4 = tmp2[[i]][tmp3]
    }
    return(tmp4)
  })
  top_linked_gene = unlist(top_linked_gene_list)
  results$top_linked_gene = top_linked_gene
  results$gene_anno = results$top_linked_gene
  results[is.na(results$gene_anno),]$gene_anno = paste0(results[is.na(results$gene_anno),]$nearest_gene, ' (nearest)')

  p = ggplot(results, aes(log2FoldChange, -log10(padj))) + geom_point(aes(col = Expression)) +
    geom_text_repel(max.overlaps = Inf, min.segment.length = -1, aes(label=c(results$gene_anno[1:20], rep('',nrow(results)-20))))+
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
  p <- vol_plot_new(peak_p2g, celltype, 'DA peaks')
  pdf(paste0('../plots/pseudo_bulk/ATAC_updated/new/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), vol_plot_new(peak_p2g, celltype, 'DA peaks'))
}

pdf(paste0('../plots/pseudo_bulk/ATAC_updated/new/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 35, height = 28)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, ncol = 5))
dev.off()

rna_list = lapply(levels[1:20], function(i){
  results <- ipa(rna, i)
  results = results[,-c('de_family','de_method','de_type')]
  return(results)
})
rna_new = Reduce(rbind, rna_list)
write.xlsx(rna_new, file = '../output/pseudo_bulk_updated/DE_rna_DESeq2_paired_new.xlsx', colNames = T, rowNames = F)

adt_list = lapply(levels[1:20], function(i){
  results <- ipa(adt, i)
  results = results[,-c('de_family','de_method','de_type')]
  return(results)
})
adt_new = Reduce(rbind, adt_list)
write.xlsx(adt_new, file = '../output/pseudo_bulk_updated/DE_adt_DESeq2_paired_new.xlsx', colNames = T, rowNames = F)

peak_anno = function(de, name){
  data <- de[de$cell_type == name,]
  results = mutate(data, Expression = ifelse(data$padj < 0.05, "DEGs", "Not Significant"))
  results = results[order(results$pvalue),]
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$log2FoldChange > 0) & (results$padj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$log2FoldChange < 0) & (results$padj < 0.05))] <- "Down-regulated"
  
  tmp1 = str_split(results$FDR, '\\|')
  index = lapply(tmp1, function(i){which.min(as.numeric(i))})
  tmp2 = str_split(results$gene, '\\|')

  top_linked_gene_list = lapply(1:length(index), function(i){
    tmp3 = index[[i]]
    if(length(tmp3) == 0){
      tmp4 = NA
    }
    else{
      tmp4 = tmp2[[i]][tmp3]
    }
    return(tmp4)
  })
  top_linked_gene = unlist(top_linked_gene_list)
  results$top_linked_gene = top_linked_gene
  results$gene_anno = results$top_linked_gene
  results[is.na(results$gene_anno),]$gene_anno = paste0(results[is.na(results$gene_anno),]$nearest_gene, ' (nearest)')
  return(results)
}
peak_list = lapply(levels[1:20], function(i){
  results <- peak_anno(peak_p2g, i)
  results = results[,c('cell_type','baseMean','log2FoldChange','pvalue', 'padj', 'peak', 'Expression', 'gene_anno', 'top_linked_gene', 'nearest_gene', 'Correlation', 'FDR', 'gene')]
  return(results)
})
peak_new = Reduce(rbind, peak_list)
write.xlsx(peak_new, file = '../output/pseudo_bulk_updated/DE_peak_DESeq2_paired_new.xlsx', colNames = T, rowNames = F)

peak_p2g = read.xlsx('../output/colocalization/peak_snp_p2g_gene_sig_consistent_updated.xlsx')
peak_p2g_consistent = peak_p2g[peak_p2g$gene_consistent_num > 0,]

ipa_p2g = function(p2g, de, name){
  data <- p2g[p2g$cell_type == name,]
  
  genes = str_split(data$gene_consistent, '\\|')
  peaks = data$peak
  pairs = lapply(1:nrow(data), function(i){
    pair = paste(peaks[i], genes[[i]], sep = ':')
    return(pair)
  })
  pairs_df = unique(unlist(pairs))
  df = data.frame(pair = pairs_df, peak = str_split_fixed(pairs_df, ':', 2)[,1], gene = str_split_fixed(pairs_df, ':', 2)[,2])

  df_list = split(df, df$gene)

  df_list_re <- as.data.frame(t(sapply(df_list, function(df){
    table <- apply(df, 2, function(i){paste0(i, collapse = '|')})
    return(table)
  })))
  df_list_re$gene <- rownames(df_list_re)
  df_list_re = df_list_re[,c('gene', 'peak')]
  df_list_re$peak_num = unlist(lapply(str_split(df_list_re$peak, '\\|'), function(i){return(length(i))}))

  de <- de[de$cell_type == name,]
  results = mutate(de, Expression = ifelse(de$padj < 0.05, "DEGs", "Not Significant"))
  results = results[order(results$pvalue),]
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$log2FoldChange > 0) & (results$padj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$log2FoldChange < 0) & (results$padj < 0.05))] <- "Down-regulated"
  
  results_sig = results[results$gene %in% genes,]
  df_list_re = df_list_re[results_sig$gene,]
  results_sig$peak_num = df_list_re$peak_num
  results_sig$peak = df_list_re$peak
  
  return(results_sig)
}

for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  results <- ipa_p2g(peak_p2g_consistent, rna, celltype)
  write.xlsx(results, file = paste0('../plots/pseudo_bulk/RNA_updated/IPA_consistent/', celltype, '.xlsx'), colNames = T, rowNames = F)
}

vol_plot_p2g <- function(p2g, de, name, title){
  data <- p2g[p2g$cell_type == name,]
  genes = unique(unlist(str_split(data$gene_consistent, '\\|')))

  data <- de[de$cell_type == name,]
  results = mutate(data, Expression = ifelse(data$padj < 0.05, "DEGs", "Not Significant"))
  results = results[order(results$pvalue),]
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$log2FoldChange > 0) & (results$padj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$log2FoldChange < 0) & (results$padj < 0.05))] <- "Down-regulated"

  results$Expression_re = "Not Significant"
  results$Expression_re[which((results$Expression == "Up-regulated") & (results$gene %in% genes))] <- "Up-regulated"
  results$Expression_re[which((results$Expression == "Down-regulated") & (results$gene %in% genes))] <- "Down-regulated"

  p = ggplot(results, aes(log2FoldChange, -log10(padj))) + geom_point(aes(col = Expression_re)) +
    geom_text_repel(max.overlaps = Inf, min.segment.length = -1, aes(label=ifelse(gene %in% results[results$Expression != "Not Significant",]$gene[1:20],gene,'')))+
    scale_color_manual(values = c("Down-regulated" = "dodgerblue3", "Not Significant" = "darkgrey", "Up-regulated" = "firebrick")) +
    theme_bw(base_size = 12) +
    theme(text = element_text(size=16)) +
    theme(legend.position = "bottom") + geom_hline(yintercept = -log10(0.05), color = "grey1", lty = 2, lwd = 0.5)+
    ylab('-log10(FDR-adjusted P-value)')+
    xlab('log2(Fold Change)')+
    labs(col = '', title = name, subtitle = paste0(ifelse(is.na(table(results$Expression)["Up-regulated"]), 0, table(results$Expression)["Up-regulated"])+ifelse(is.na(table(results$Expression)["Down-regulated"]), 0, table(results$Expression)["Down-regulated"]), ' significant ', title, ' were identified'))+
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
  return(p)
}

for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  p <- vol_plot_p2g(peak_p2g_consistent, rna, celltype, 'DE genes')
  pdf(paste0('../plots/pseudo_bulk/RNA_updated/new/', celltype, '.pdf'), width = 7, height = 7)
  print(p)
  dev.off()
  assign(paste0('p',i), vol_plot_p2g(peak_p2g_consistent, rna, celltype, 'DE genes'))
}

pdf(paste0('../plots/pseudo_bulk/RNA_updated/new/', 'Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23', '.pdf'), width = 35, height = 28)
print(wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, ncol = 5))
dev.off()

####################
rna_list = lapply(levels[1:20], function(i){
  results <- ipa(rna, i)
  results = results[,-c('de_family','de_method','de_type')]
  return(results)
})
rna_new = Reduce(rbind, rna_list)
write.xlsx(rna_new, file = '../output/pseudo_bulk_updated/DE_rna_DESeq2_paired_new.xlsx', colNames = T, rowNames = F)

adt_list = lapply(levels[1:20], function(i){
  results <- ipa(adt, i)
  results = results[,-c('de_family','de_method','de_type')]
  return(results)
})
adt_new = Reduce(rbind, adt_list)
write.xlsx(adt_new, file = '../output/pseudo_bulk_updated/DE_adt_DESeq2_paired_new.xlsx', colNames = T, rowNames = F)

peak_anno = function(de, name){
  data <- de[de$cell_type == name,]
  results = mutate(data, Expression = ifelse(data$padj < 0.05, "DEGs", "Not Significant"))
  results = results[order(results$pvalue),]
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$log2FoldChange > 0) & (results$padj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$log2FoldChange < 0) & (results$padj < 0.05))] <- "Down-regulated"
  
  tmp1 = str_split(results$FDR, '\\|')
  index = lapply(tmp1, function(i){which.min(as.numeric(i))})
  tmp2 = str_split(results$gene, '\\|')
  
  top_linked_gene_list = lapply(1:length(index), function(i){
    tmp3 = index[[i]]
    if(length(tmp3) == 0){
      tmp4 = NA
    }
    else{
      tmp4 = tmp2[[i]][tmp3]
    }
    return(tmp4)
  })
  top_linked_gene = unlist(top_linked_gene_list)
  results$top_linked_gene = top_linked_gene
  results$gene_anno = results$top_linked_gene
  results[is.na(results$gene_anno),]$gene_anno = paste0(results[is.na(results$gene_anno),]$nearest_gene, ' (nearest)')
  return(results)
}
peak_list = lapply(levels[1:20], function(i){
  results <- peak_anno(peak_p2g, i)
  results = results[,c('cell_type','baseMean','log2FoldChange','pvalue', 'padj', 'peak', 'Expression', 'gene_anno', 'top_linked_gene', 'nearest_gene', 'Correlation', 'FDR', 'gene')]
  return(results)
})
peak_new = Reduce(rbind, peak_list)
write.xlsx(peak_new, file = '../output/pseudo_bulk_updated/DE_peak_DESeq2_paired_new.xlsx', colNames = T, rowNames = F)

####################
rna_list = lapply(levels[1:20], function(i){
  results <- ipa(rna, i)
  results = table(results$Expression)[c('Up-regulated', 'Down-regulated', 'Not Significant')]
  return(results)
})
rna_new = as.data.frame(Reduce(rbind, rna_list))
rownames(rna_new) <- levels
rna_new$Sum <- rna_new$`Down-regulated`+rna_new$`Up-regulated`
data_plot <- rna_new[,c('Up-regulated', 'Down-regulated')]
data_plot$id <- rownames(data_plot)
data_plot <- reshape2::melt(data_plot)

library(scales)
p1 = ggplot(data_plot, aes(x = value, y = id, fill = variable))+
  geom_col(position = 'stack', col = 'white')+
  geom_text(aes(label = id, x = 1e3), hjust = 0, colour = "black")+
  theme_bw()+
  scale_fill_manual(values = c('Up-regulated' = 'firebrick', 'Down-regulated' = 'dodgerblue3'))+
  scale_x_continuous(n.breaks = 12)+
  labs(y = 'Cell types', x = '# of DE genes', fill = '')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'bottom', axis.text.y = element_blank(), axis.ticks.y = element_blank())
pdf('../output/pseudo_bulk_updated/num_DE_genes.pdf', width = 3, height = 5)
p1
dev.off()


peak_list = lapply(levels[1:20], function(i){
  results <- ipa(peak, i)
  results = table(results$Expression)[c('Up-regulated', 'Down-regulated', 'Not Significant')]
  results[is.na(results)] <- 0
  return(results)
})
peak_new = as.data.frame(Reduce(rbind, peak_list))
rownames(peak_new) <- levels
peak_new$Sum <- peak_new$`Down-regulated`+peak_new$`Up-regulated`
data_plot <- peak_new[,c('Up-regulated', 'Down-regulated')]
data_plot$id <- rownames(data_plot)
data_plot <- reshape2::melt(data_plot)

p1 = ggplot(data_plot, aes(x = value, y = id, fill = variable))+
  geom_col(position = 'stack', col = 'white')+
  geom_text(aes(label = id, x = 1e3), hjust = 0, colour = "black")+
  theme_bw()+
  scale_fill_manual(values = c('Up-regulated' = 'firebrick', 'Down-regulated' = 'dodgerblue3'))+
  scale_x_continuous(n.breaks = 12)+
  labs(y = 'Cell types', x = '# of DA ATAC peaks', fill = '')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'bottom', axis.text.y = element_blank(), axis.ticks.y = element_blank())
pdf('../output/pseudo_bulk_updated/num_DA_peaks.pdf', width = 3, height = 5)
p1
dev.off()
