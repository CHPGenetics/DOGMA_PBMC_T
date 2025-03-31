library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
library(ggpubr)
library(cowplot)
library(patchwork)
library(DESeq2)
library(dplyr)
library(ggplot2)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')

df = read.xlsx('../output/colocalization/peak_snp_p2g_gene_sig_consistent_only_snp_updated.xlsx')
table(unlist(str_split(df$gene_consistent, '\\|')))[order(table(unlist(str_split(df$gene_consistent, '\\|'))), decreasing = T)[1:20]]
# PTGER4      TTC33  SATB1-AS1 AC144521.1      SATB1       CREM     GPR183      LNPEP        MYC    RASGRP1 
# 17         14         13         10         10          7          7          7          7          7 
# TNFRSF14      IL2RA  LINC00824     MAP3K8      TAGAP AC009126.1       CAST       ELL2        LIF      RBM17 
# 7          6          6          6          6          5          5          5          5          5 

check_gene <- function(df, gene){
  genes <- str_split(df$gene_consistent, '\\|')
  tmp <- sapply(genes, function(i){any(i %in% gene)})
  return(tmp)
}

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

draw_plot <- function(peak, gene, celltype, snp){
  atac_bulk <- readRDS(paste0('../output/eQTA_updated/matrix_2/ATAC/', celltype, '.RDS'))
  rna_bulk <- readRDS(paste0('../output/eQTA_updated/matrix_2/RNA/', celltype, '.RDS'))
  
  shared = intersect(colnames(atac_bulk), colnames(rna_bulk))

  atac_bulk = atac_bulk[,shared]
  rna_bulk = rna_bulk[,shared]
  condition = factor(str_split_fixed(shared, ':', 2)[,2], levels = c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'))
  sample = str_split_fixed(shared, ':', 2)[,1]
  
  data_plot <- data.frame(rna = rna_bulk[gene,], atac = atac_bulk[peak,], condition = condition, sample = sample)
  
  check = names(which(table(data_plot$sample) == 1))
  if (length(check) > 0){
    data_plot = data_plot[!data_plot$sample %in% check,]
  }

  p1 = ggscatter(data_plot, x = 'atac', y = 'rna', color = 'condition',
                  rug = T,
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.x = min(data_plot$atac), label.y = max(data_plot$rna)*1.1 - min(data_plot$rna)*0.1)) +
    scale_color_manual(values = c('firebrick', 'dodgerblue3'))+
    labs(x = peak, y = gene)+
    labs(title = celltype, col = '')+
    theme_cowplot()+
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom')
  
  p2 = ggpaired(data_plot, x = 'condition', y = 'rna', fill = 'condition', id = 'sample',
    line.color = "darkgray", line.size = 0.4)+
    scale_fill_manual(values = c('firebrick', 'dodgerblue3'))+
    stat_compare_means(method = 'wilcox.test', label = "p.signif", paired = TRUE, comparisons = list(c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2")))+
    labs(x = '', y = gene)+
    labs(fill = '')+
    theme_cowplot()+
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none', axis.text.x = element_blank())+ylim(min(data_plot$rna), max(data_plot$rna)*1.1 - min(data_plot$rna)*0.1)

  p3 = ggpaired(data_plot, x = 'condition', y = 'atac', fill = 'condition', id = 'sample',
    line.color = "darkgray", line.size = 0.4)+
    scale_fill_manual(values = c('firebrick', 'dodgerblue3'))+
    stat_compare_means(method = 'wilcox.test', label = "p.signif", paired = TRUE, comparisons = list(c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2")))+
    labs(x = '', y = peak)+
    labs(title = paste(snp, collapse = '\n'), fill = '')+
    theme_cowplot()+
    theme(plot.title = element_text(size = 10, hjust = 0.5), legend.position = 'none', axis.text.x = element_blank())+ylim(min(data_plot$atac), max(data_plot$atac)*1.1 - min(data_plot$atac)*0.1)

  p = p1 | (p2|p3)
  return(p)
}

View(df[check_gene(df, 'PTGER4'),])
df0 = df[check_gene(df, 'PTGER4'),]

for(i in 1:nrow(df0)){
    gene_name = 'PTGER4'
    peak_name = df0$peak[i]
    cell_type = df0$cell_type[i]
    snp = str_split(df0$SNP[i], '\\|')[[1]]
    p = draw_plot(peak_name, gene_name, cell_type, snp)
    pdf(paste0('../plots/visualization_2/', paste(gene_name,peak_name,cell_type,sep = '-'), '.pdf'), width = 7.5, height = 5)
    print(p)
    dev.off()
}

View(df[check_gene(df, 'TTC33'),])
df0 = df[check_gene(df, 'TTC33'),]

for(i in 1:nrow(df0)){
    gene_name = 'TTC33'
    peak_name = df0$peak[i]
    cell_type = df0$cell_type[i]
    snp = str_split(df0$SNP[i], '\\|')[[1]]
    p = draw_plot(peak_name, gene_name, cell_type, snp)
    pdf(paste0('../plots/visualization_2/', paste(gene_name,peak_name,cell_type,sep = '-'), '.pdf'), width = 7.5, height = 5)
    print(p)
    dev.off()
}

View(df[check_gene(df, 'SATB1-AS1'),])
df0 = df[check_gene(df, 'SATB1-AS1'),]

for(i in 1:nrow(df0)){
    gene_name = 'SATB1-AS1'
    peak_name = df0$peak[i]
    cell_type = df0$cell_type[i]
    snp = str_split(df0$SNP[i], '\\|')[[1]]
    p = draw_plot(peak_name, gene_name, cell_type, snp)
    pdf(paste0('../plots/visualization_2/', paste(gene_name,peak_name,cell_type,sep = '-'), '.pdf'), width = 7.5, height = 5)
    print(p)
    dev.off()
}

View(df[check_gene(df, 'SATB1'),])
df0 = df[check_gene(df, 'SATB1'),]

for(i in 1:nrow(df0)){
    gene_name = 'SATB1'
    peak_name = df0$peak[i]
    cell_type = df0$cell_type[i]
    snp = str_split(df0$SNP[i], '\\|')[[1]]
    p = draw_plot(peak_name, gene_name, cell_type, snp)
    pdf(paste0('../plots/visualization_2/', paste(gene_name,peak_name,cell_type,sep = '-'), '.pdf'), width = 7.5, height = 5)
    print(p)
    dev.off()
}

View(df[check_gene(df, 'CREM'),])
df0 = df[check_gene(df, 'CREM'),]

for(i in 1:nrow(df0)){
    gene_name = 'CREM'
    peak_name = df0$peak[i]
    cell_type = df0$cell_type[i]
    snp = str_split(df0$SNP[i], '\\|')[[1]]
    p = draw_plot(peak_name, gene_name, cell_type, snp)
    pdf(paste0('../plots/visualization_2/', paste(gene_name,peak_name,cell_type,sep = '-'), '.pdf'), width = 7.5, height = 5)
    print(p)
    dev.off()
}

View(df[check_gene(df, 'GPR183'),])
df0 = df[check_gene(df, 'GPR183'),]

for(i in 1:nrow(df0)){
  gene_name = 'GPR183'
  peak_name = df0$peak[i]
  cell_type = df0$cell_type[i]
  snp = str_split(df0$SNP[i], '\\|')[[1]]
  p = draw_plot(peak_name, gene_name, cell_type, snp)
  pdf(paste0('../plots/visualization_2/', paste(gene_name,peak_name,cell_type,sep = '-'), '.pdf'), width = 7.5, height = 5)
  print(p)
  dev.off()
}

View(df[check_gene(df, 'LNPEP'),])
df0 = df[check_gene(df, 'LNPEP'),]

for(i in 1:nrow(df0)){
    gene_name = 'LNPEP'
    peak_name = df0$peak[i]
    cell_type = df0$cell_type[i]
        snp = str_split(df0$SNP[i], '\\|')[[1]]
    p = draw_plot(peak_name, gene_name, cell_type, snp)
    pdf(paste0('../plots/visualization_2/', paste(gene_name,peak_name,cell_type,sep = '-'), '.pdf'), width = 7.5, height = 5)
    print(p)
    dev.off()
}

View(df[check_gene(df, 'MYC'),])
df0 = df[check_gene(df, 'MYC'),]

for(i in 1:nrow(df0)){
  gene_name = 'MYC'
  peak_name = df0$peak[i]
  cell_type = df0$cell_type[i]
  snp = str_split(df0$SNP[i], '\\|')[[1]]
  p = draw_plot(peak_name, gene_name, cell_type, snp)
  pdf(paste0('../plots/visualization_2/', paste(gene_name,peak_name,cell_type,sep = '-'), '.pdf'), width = 7.5, height = 5)
  print(p)
  dev.off()
}

View(df[check_gene(df, 'RASGRP1'),])
df0 = df[check_gene(df, 'RASGRP1'),]

for(i in 1:nrow(df0)){
  gene_name = 'RASGRP1'
  peak_name = df0$peak[i]
  cell_type = df0$cell_type[i]
  snp = str_split(df0$SNP[i], '\\|')[[1]]
  p = draw_plot(peak_name, gene_name, cell_type, snp)
  pdf(paste0('../plots/visualization_2/', paste(gene_name,peak_name,cell_type,sep = '-'), '.pdf'), width = 7.5, height = 5)
  print(p)
  dev.off()
}

View(df[check_gene(df, 'TNFRSF14'),])
df0 = df[check_gene(df, 'TNFRSF14'),]

for(i in 1:nrow(df0)){
    gene_name = 'TNFRSF14'
    peak_name = df0$peak[i]
    cell_type = df0$cell_type[i]
        snp = str_split(df0$SNP[i], '\\|')[[1]]
    p = draw_plot(peak_name, gene_name, cell_type, snp)
    pdf(paste0('../plots/visualization_2/', paste(gene_name,peak_name,cell_type,sep = '-'), '.pdf'), width = 7.5, height = 5)
    print(p)
    dev.off()
}

View(df[check_gene(df, 'IL2RA'),])
df0 = df[check_gene(df, 'IL2RA'),]

for(i in 1:nrow(df0)){
    gene_name = 'IL2RA'
    peak_name = df0$peak[i]
    cell_type = df0$cell_type[i]
        snp = str_split(df0$SNP[i], '\\|')[[1]]
    p = draw_plot(peak_name, gene_name, cell_type, snp)
    pdf(paste0('../plots/visualization_2/', paste(gene_name,peak_name,cell_type,sep = '-'), '.pdf'), width = 7.5, height = 5)
    print(p)
    dev.off()
}

