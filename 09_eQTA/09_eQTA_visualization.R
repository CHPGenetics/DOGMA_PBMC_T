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
library(stringr)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')

p2g_bp <- read.xlsx('../output/colocalization/p2g_snp/p2g_snp_immune_sub.xlsx')

load('../output/ArchR/DOGMA_filtered_multiome/PeakMatrix_grange.RData')
peak_df <- as.data.frame(peak_gr)
peak_df$peak <- paste(peak_df$seqnames, peak_df$start, peak_df$end, sep = '-')

draw_plot <- function(peak, gene, snp, celltype){
  # atac <- readRDS(paste0('../output/ArchR/DOGMA_filtered_multiome/Peak2GeneLinks/', celltype, '/seATAC-Group-KNN.rds'))
  # rna <- readRDS(paste0('../output/ArchR/DOGMA_filtered_multiome/Peak2GeneLinks/', celltype, '/seRNA-Group-KNN.rds'))
  # 
  # atac_mat <- assays(atac)[[1]]
  # rownames(atac_mat) <- peak_df$peak
  # colnames(atac_mat) <- paste0('meta', 1:dim(atac_mat)[2])
  # 
  # rna_mat <- assays(rna)[[1]]
  # rownames(rna_mat) <- rowData(rna)$name
  # colnames(rna_mat) <- paste0('meta', 1:dim(rna_mat)[2])
  
  atac_bulk <- readRDS(paste0('../output/eQTA/matrix/ATAC/', celltype, '.RDS'))
  rna_bulk <- readRDS(paste0('../output/eQTA/matrix/RNA/', celltype, '.RDS'))
  
  # data_plot <- data.frame(rna = rna_mat[gene,], atac = atac_mat[peak,])
  # data_plot <- data_plot[data_plot$atac > 0,]
  # p1 <- ggscatter(data_plot, x = 'atac', y = 'rna',
  #                 add = "reg.line",  # Add regressin line
  #                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
  #                 conf.int = TRUE, # Add confidence interval
  #                 cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
  #                 cor.coeff.args = list(method = "pearson", label.x = min(data_plot$atac), label.y = max(data_plot$rna)*1.2)) +
  #   labs(x = peak, y = gene)+
  #   ggtitle('Metacell')+
  #   theme_cowplot()+
  #   theme(plot.title = element_text(hjust = 0.5))
  
  data_plot <- data.frame(rna = rna_bulk[gene,], atac = atac_bulk[peak,])
  p2 <- ggscatter(data_plot, x = 'atac', y = 'rna',
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.x = min(data_plot$atac), label.y = max(data_plot$rna)*1.2)) +
    labs(x = peak, y = gene)+
    labs(title = celltype)+
    theme_cowplot()+
    theme(plot.title = element_text(hjust = 0.5))
  
  p <- wrap_plots(p2, grid::textGrob(snp), ncol = 2)
  return(p)
}

filter_p2g <- function(p2g, disease, gene){
  tmp <- p2g[str_detect(p2g$Disease_merged, disease),]
  tmp <- tmp[tmp$gene == gene,]
  return(tmp)
}

p2g <- filter_p2g(p2g_bp, 'Psoriasis', 'IFNLR1')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/Psoriasis_IFNLR1.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'Psoriasis', 'RUNX3')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/Psoriasis_RUNX3.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'Behcets_disease', 'KLRC4')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/Behcets_disease_KLRC4.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'Alopecia_areata', 'IL21')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/Alopecia_areata_IL21.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'Ankylosing_spondylitis', 'GPR65')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/Ankylosing_spondylitis_GPR65.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'Systemic_lupus_erythematosus', 'BLK')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/Systemic_lupus_erythematosus_BLK.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'Systemic_lupus_erythematosus', 'CD44')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/Systemic_lupus_erythematosus_CD44.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g_bp <- read.xlsx('../output/colocalization/p2g_snp/p2g_snp_ibd_sub.xlsx')
p2g_bp <- p2g_bp[!is.na(p2g_bp$SNP),]

p2g <- filter_p2g(p2g_bp, 'All', 'PTGER4')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/IBD_PTGER4.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'All', 'LNPEP')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/IBD_LNPEP.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'All', 'TNFRSF14')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/IBD_TNFRSF14.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'All', 'IL2RA')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/IBD_IL2RA.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'All', 'CCL20')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/IBD_CCL20.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'All', 'CXCR5')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/IBD_CXCR5.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'All', 'ITGAL')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/IBD_ITGAL.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############
##############

p2g <- filter_p2g(p2g_bp, 'Lead', 'CCL20')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/IBD_CCL20_lead.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'Lead', 'CXCR5')

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/IBD_CXCR5_lead.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############

p2g <- filter_p2g(p2g_bp, 'Lead', 'ITGAL')
p2g$SNP <- 'rs11150589'

p_list <- list()
for(i in 1:nrow(p2g)){
  peak <- p2g$peak[i]
  gene <- p2g$gene[i]
  snp <- paste0(str_split(p2g$SNP[i], '\\|')[[1]], collapse = '\n')
  celltype <- p2g$celltype[i]
  
  p_list[[i]] <- draw_plot(peak, gene, snp, celltype)
}

pdf('../plots/eQTA_replicate_0.05/SNP/IBD_ITGAL_lead.pdf', height = ifelse(length(p_list) %% 4 == 0, 4*(length(p_list) %/% 4), 4*(length(p_list) %/% 4 + 1)), width = 32)
wrap_plots(p_list, ncol = 4)
dev.off()
##############