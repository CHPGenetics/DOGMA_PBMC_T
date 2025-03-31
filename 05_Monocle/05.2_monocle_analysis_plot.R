library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(data.table)
library(Matrix)
library(dplyr)
library(stringr)
library(future)
library(harmony)
library(clustree)
library(openxlsx)
library(monocle3)
library(SeuratWrappers)
library(magrittr)
library(ComplexHeatmap)
library(ggsci)
library(cowplot)
library(splines)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

data <- readRDS('../output/tcell_annotated_updated.RDS')
add_pseudotime <- function(name){
  pseudotime <- read.csv(paste0('../output/monocle/', name, '.csv'), row.names = 'X')
  data <- AddMetaData(data, pseudotime, paste0('monocle_', name))
  return(data)
}
data <- add_pseudotime('cd4_naive')
data <- add_pseudotime('cd4_memory')
data <- add_pseudotime('cd8_naive')
data <- add_pseudotime('cd8_memory')
data <- add_pseudotime('mait')
# data <- add_pseudotime('Th17_Th1')

motif <- readRDS('../output/ArchR/DOGMA_filtered_multiome/motif.RDS')
motif <- assays(motif)[[2]]
colnames(motif) <- gsub('#', '_', colnames(motif))
motif <- motif[,colnames(data)]
rownames(motif) <- gsub('_', '-', rownames(motif))
data[['chromvar2']] <- CreateAssayObject(data = motif)
rm(motif)

rna <- t(data@assays$SCT@data)
adt <- t(data@assays$ADT@data)
motif <- t(data@assays$chromvar2@data)
meta <- data@meta.data

marker <- read.csv('../output/ArchR/DOGMA_filtered_multiome/Plots/Plot-MonocleCD4Naive_fine-Traj-Heatmap_Cor.csv')

plot_marker <- function(type, marker){
  rna0 <- rna[,marker]
  motif0 <- motif[,colnames(motif)[grep(paste0('^', marker, '-', collapse = '|'), colnames(motif))]]
  mtx0 <- meta[,c(type, 'celltype', 'condition')]
  cell <- mtx0[!is.na(mtx0[,1]),]
  
  rna0 <- rna0[rownames(cell),]
  # rna0$mean <- rowMeans(rna0)
  rna0_s <- scale(as.matrix(rna0))
  # rna0_s$mean <- rowMeans(rna0_s)
  motif0 <- motif0[rownames(cell),]
  # motif0$mean <- rowMeans(motif0)
  motif0_s <- scale(as.matrix(motif0))
  # motif0_s$mean <- rowMeans(motif0_s)
  mtx0 <- mtx0[rownames(cell),]
  
  tmp1 <- cbind(mtx0, rna0_s)
  tmp1 <- melt(tmp1, id.vars = c(type, 'celltype', 'condition'))
  tmp1$variable <- factor(tmp1$variable, levels = c(marker))
  
  tmp2 <- cbind(mtx0, motif0_s)
  tmp2 <- melt(tmp2, id.vars = c(type, 'celltype', 'condition'))
  tmp2$variable_order <- str_split_fixed(tmp2$variable, '-', 2)[,1]
  tmp2$variable_order <- factor(tmp2$variable_order, levels = c(marker))
  # ggplot(cell, aes(x = monocle_cd4_naive)) +
  #   geom_density(alpha = .4, aes(fill = celltype), col = "transparent") +
  #   geom_density(aes(col = celltype), fill = "transparent", size = 1.5) +
  #   guides(col = FALSE) +
  #   scale_fill_npg() +
  #   scale_color_npg() +
  #   labs(x = "Pseudotime", y = "Density", fill = "")+
  #   theme_cowplot()
  
  p1 <- ggplot(tmp1, aes(x = get(type), y = value, col = variable)) +
    geom_smooth(formula = y ~ s(x, bs = "cs"), se = T, alpha = 0.2)+
    geom_smooth(formula = y ~ s(x, bs = "cs"), se = F, alpha = 0.2, aes(x = get(type), y = value), linetype = 'dashed', col = 'black')+
    scale_color_npg() +
    labs(x = "Pseudotime", y = "", col = "")+
    theme_cowplot()+
    theme(legend.position = 'none')
  
  p2 <- ggplot(tmp2, aes(x = get(type), y = value, col = variable_order)) +
    geom_smooth(formula = y ~ s(x, bs = "cs"), se = T, alpha = 0.2)+
    geom_smooth(formula = y ~ s(x, bs = "cs"), se = F, alpha = 0.2, aes(x = get(type), y = value), linetype = 'dashed', col = 'black')+
    scale_color_npg() +
    labs(x = "Pseudotime", y = "", col = "")+
    theme_cowplot()
  
  plot <- list(p1, p2)
  return(plot)
}

plot_pseudotime <- function(type, density){
  # rna0 <- rna[,marker]
  # motif0 <- motif[,colnames(motif)[grep(paste0('^', marker, collapse = '|'), colnames(motif))]]
  mtx0 <- meta[,c(type, 'celltype', 'condition')]
  cell <- mtx0[!is.na(mtx0[,1]),]
  
  # rna0 <- as.data.frame(as.matrix(rna0[rownames(cell),]))
  # rna0$mean <- rowMeans(rna0)
  # rna0_s <- scale(as.matrix(rna0))
  # motif0 <- as.data.frame(as.matrix(motif0[rownames(cell),]))
  # motif0$mean <- rowMeans(motif0)
  # motif0_s <- scale(as.matrix(motif0))
  mtx0 <- mtx0[rownames(cell),]
  
  # tmp1 <- cbind(mtx0, rna0_s)
  # tmp1 <- melt(tmp1, id.vars = c(type, 'celltype', 'condition'))
  # tmp1$variable <- factor(tmp1$variable, levels = c(marker, 'mean'))
  # 
  # tmp2 <- cbind(mtx0, motif0_s)
  # tmp2 <- melt(tmp2, id.vars = c(type, 'celltype', 'condition'))
  # tmp2$variable_order <- str_split_fixed(tmp2$variable, '-', 2)[,1]
  # tmp2$variable_order <- factor(tmp2$variable_order, levels = c(marker, 'mean'))
  p <- ggplot(cell, aes(x = get(type))) +
    geom_density(alpha = .4, aes(fill = get(density)), col = "transparent") +
    geom_density(aes(col = get(density)), fill = "transparent", size = 1.5) +
    guides(col = FALSE) +
    scale_fill_npg() +
    scale_color_npg() +
    labs(x = "Pseudotime", y = "Density", fill = "")+
    theme_cowplot()+
    theme(legend.position = 'bottom')
  
  # p1 <- ggplot(tmp1[tmp1$variable != 'mean',], aes(x = get(type), y = value, col = variable)) +
  #   geom_smooth(formula = y ~ ns(x, df = 3), se = T, alpha = 0.2)+
  #   geom_smooth(formula = y ~ ns(x, df = 3), se = F, alpha = 0.2, data = tmp1[tmp1$variable == 'mean',], linetype = 'dashed', col = 'black')+
  #   scale_color_npg() +
  #   labs(x = "Pseudotime", y = "", col = "")+
  #   theme_cowplot()+
  #   theme(legend.position = 'none')
  # 
  # p2 <- ggplot(tmp2[tmp2$variable_order != 'mean',], aes(x = get(type), y = value, col = variable_order)) +
  #   geom_smooth(formula = y ~ ns(x, df = 3), se = T, alpha = 0.2)+
  #   geom_smooth(formula = y ~ ns(x, df = 3), se = F, alpha = 0.2, data = tmp2[tmp2$variable_order == 'mean',], linetype = 'dashed', col = 'black')+
  #   scale_color_npg() +
  #   labs(x = "Pseudotime", y = "", col = "")+
  #   theme_cowplot()
  
  # plot <- list(p1, p2)
  return(p)
}

plot_coexpresion <- function(type, marker){
  rna0 <- rna[,marker]
  motif0 <- motif[,colnames(motif)[grep(paste0('^', marker, '-', collapse = '|'), colnames(motif))]]
  mtx0 <- meta[,c(type, 'celltype', 'condition')]
  cell <- mtx0[!is.na(mtx0[,1]),]
  
  rna0 <- as.matrix(rna0[rownames(cell),])
  # rna0$mean <- rowMeans(rna0)
  # rna0_s <- scale(as.matrix(rna0))
  motif0 <- as.matrix(motif0[rownames(cell),])
  # motif0$mean <- rowMeans(motif0)
  # motif0_s <- scale(as.matrix(motif0))
  mtx0 <- mtx0[rownames(cell),]
  
  tmp1 <- cor(rna0)
  tmp2 <- cor(motif0)
  # tmp1 <- cbind(mtx0, rna0_s)
  # tmp1 <- melt(tmp1, id.vars = c(type, 'celltype', 'condition'))
  # tmp1$variable <- factor(tmp1$variable, levels = c(marker, 'mean'))
  # 
  # tmp2 <- cbind(mtx0, motif0_s)
  # tmp2 <- melt(tmp2, id.vars = c(type, 'celltype', 'condition'))
  # tmp2$variable_order <- str_split_fixed(tmp2$variable, '-', 2)[,1]
  # tmp2$variable_order <- factor(tmp2$variable_order, levels = c(marker, 'mean'))
  p1 <- pheatmap(tmp1, cluster_rows = T, cluster_cols = T, show_colnames = F)
  p2 <- pheatmap(tmp2, cluster_rows = T, cluster_cols = T, show_colnames = F)
  
  # p1 <- ggplot(tmp1[tmp1$variable != 'mean',], aes(x = get(type), y = value, col = variable)) +
  #   geom_smooth(formula = y ~ ns(x, df = 3), se = T, alpha = 0.2)+
  #   geom_smooth(formula = y ~ ns(x, df = 3), se = F, alpha = 0.2, data = tmp1[tmp1$variable == 'mean',], linetype = 'dashed', col = 'black')+
  #   scale_color_npg() +
  #   labs(x = "Pseudotime", y = "", col = "")+
  #   theme_cowplot()+
  #   theme(legend.position = 'none')
  # 
  # p2 <- ggplot(tmp2[tmp2$variable_order != 'mean',], aes(x = get(type), y = value, col = variable_order)) +
  #   geom_smooth(formula = y ~ ns(x, df = 3), se = T, alpha = 0.2)+
  #   geom_smooth(formula = y ~ ns(x, df = 3), se = F, alpha = 0.2, data = tmp2[tmp2$variable_order == 'mean',], linetype = 'dashed', col = 'black')+
  #   scale_color_npg() +
  #   labs(x = "Pseudotime", y = "", col = "")+
  #   theme_cowplot()
  
  plot <- list(p1, p2)
  return(plot)
}

p <- plot_pseudotime('monocle_cd4_naive', 'celltype')
pdf('../plots/monocle/cd4_naive_density.pdf', width = 6, height = 6)
p
dev.off()
p <- plot_pseudotime('monocle_cd8_naive', 'celltype')
pdf('../plots/monocle/cd8_naive_density.pdf', width = 6, height = 6)
p
dev.off()
p <- plot_pseudotime('monocle_cd8_memory', 'celltype')
pdf('../plots/monocle/cd8_memory_density.pdf', width = 6, height = 6)
p
dev.off()
p <- plot_pseudotime('monocle_mait', 'celltype')
pdf('../plots/monocle/mait_density.pdf', width = 6, height = 6)
p
dev.off()

p_cor <- plot_coexpresion('monocle_cd4_naive', marker$x)
pdf('../plots/monocle/cd4_naive_cor1.pdf', width = 16, height = 14)
p_cor[[1]]
dev.off()
pdf('../plots/monocle/cd4_naive_cor2.pdf', width = 16, height = 14)
p_cor[[2]]
dev.off()

p1 <- plot_marker('monocle_cd4_naive', marker$x[marker$x %in% c("ELF2", "ETV6", "FLI1", "ETS1", "FOXO3", "FOXP1", "IKZF1", "ELF1")])
p2 <- plot_marker('monocle_cd4_naive', marker$x[marker$x %in% c("KLF5", "KLF13", "SP4", "SP3", "NFYC", "KLF2", "KLF3", "DBP", "USF2")])
p3 <- plot_marker('monocle_cd4_naive', marker$x[marker$x %in% c("MAFK", "BATF", "BATF3", "FOS", "FOSL1", "JDP2", "JUND", "BACH2", "BACH1", "NFE2L1")])
p4 <- plot_marker('monocle_cd4_naive', marker$x[marker$x %in% c("NFKB1", "NFKB2", "REL", "RELB")])
p5 <- plot_marker('monocle_cd4_naive', marker$x[marker$x %in% c("ATF4", "JUN", "CREB1")])
p6 <- plot_marker('monocle_cd4_naive', marker$x[marker$x %in% c("TBX21", "MGA")])
p7 <- plot_marker('monocle_cd4_naive', marker$x[marker$x %in% c("PRDM4", "IRF1", "IRF9", "IRF8", "IRF4")])
p8 <- plot_marker('monocle_cd4_naive', marker$x[marker$x %in% c("LEF1", "TCF7")])

pdf('../plots/monocle/cd4_naive_module.pdf', width = 8, height = 16)
wrap_plots(p1[[1]],p1[[2]],p2[[1]],p2[[2]],p3[[1]],p3[[2]],p4[[1]],p4[[2]],p5[[1]],p5[[2]],p6[[1]],p6[[2]],p7[[1]],p7[[2]],p8[[1]],p8[[2]],
           ncol = 2)
dev.off()

pdf('../plots/monocle/cd4_naive_module_updated.pdf', width = 6, height = 12)
wrap_plots(p1[[1]],p1[[2]],p2[[1]],p2[[2]],p3[[1]],p3[[2]],p4[[1]],p4[[2]],
           ncol = 2)
dev.off()
