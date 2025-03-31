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
library(cowplot)
library(ggsci)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

data <- readRDS('../output/monocle_2/tcell_annotated_updated_2_conditions_cd4_memory.RDS')
add_pseudotime <- function(name){
  pseudotime <- read.csv(paste0('../output/monocle_2/cd4_memory/', name, '.csv'), row.names = 'X')
  data <- AddMetaData(data, pseudotime, paste0('monocle_', name))
  return(data)
}

data <- add_pseudotime('Th17_Th1_1')
data <- add_pseudotime('Th17_Th1_2')
data <- add_pseudotime('Th17_Th1_6')
data_bp <- data

plot_marker_condition <- function(type, marker_type, marker, condition){
  if(marker_type == 'rna'){
    rna0 <- as.data.frame(rna[,marker])
    colnames(rna0) <- marker
    title <- paste0(marker, ' (RNA)')
  }
  else if(marker_type == 'adt'){
    rna0 <- as.data.frame(adt[,marker])
    colnames(rna0) <- marker
    title <- paste0(marker, ' (ADT)')
  }
  else if(marker_type == 'motif'){
    rna0 <- as.data.frame(motif[,colnames(motif)[grep(paste0('^', marker, '-', collapse = '|'), colnames(motif))]])
    colnames(rna0) <- marker
    title <- paste0(marker, ' (ATAC)')
  }
  
  mtx0 <- meta[,c(type, 'celltype', 'condition')]
  cell <- mtx0[!is.na(mtx0[,1]),]
  
  rna0 <- as.data.frame(rna0[rownames(cell),])
  colnames(rna0) <- marker

  mtx0 <- mtx0[rownames(cell),]
  
  tmp1 <- cbind(mtx0, rna0)
  tmp1 <- melt(tmp1, id.vars = c(type, 'celltype', 'condition'))
  
  p <- ggplot(tmp1[tmp1$condition %in% condition,], aes(x = get(type), y = value, col = condition)) +
    geom_point(size = 0.1)+
    geom_smooth(formula = y ~ s(x, bs = "cs"), se = T, alpha = 0.2)+
    geom_smooth(formula = y ~ s(x, bs = "cs"), se = F, alpha = 0.2, aes(x = get(type), y = value), linetype = 'dashed', col = 'black')+
    scale_color_npg() +
    labs(x = "Pseudotime", y = "", col = "")+
    theme_cowplot()+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom')

  return(p)
}

data <- data_bp[,!is.na(data_bp$monocle_Th17_Th1_1)]
rna <- t(data@assays$SCT@data)
adt <- t(data@assays$ADT@data)
motif <- t(data@assays$chromvar2@data)
meta <- data@meta.data

p1 <- plot_marker_condition('monocle_Th17_Th1_1', 'rna', 'IL17F', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p2 <- plot_marker_condition('monocle_Th17_Th1_1', 'rna', 'IL23R', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p3 <- plot_marker_condition('monocle_Th17_Th1_1', 'rna', 'IFNG', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p4 <- plot_marker_condition('monocle_Th17_Th1_1', 'motif', 'RORC', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p5 <- plot_marker_condition('monocle_Th17_Th1_1', 'motif', 'RORA', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'))
p6 <- plot_marker_condition('monocle_Th17_Th1_1', 'motif', 'TBX21', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p7 <- plot_marker_condition('monocle_Th17_Th1_1', 'adt', 'CD196', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p8 <- plot_marker_condition('monocle_Th17_Th1_1', 'adt', 'CD194', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p9 <- plot_marker_condition('monocle_Th17_Th1_1', 'adt', 'CD183', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_1/pathogenic_Th17.pdf', width = 12, height = 12)
wrap_plots(p1, p2, p3, p7, p8, p9, p4, p5, p6, ncol = 3)
dev.off()

p1 <- FeatureScatter(data, 'chromvar2_RORC-351', 'adt_CD196')
p1_data <- p1$data

p1_new <- ggplot(p1_data, aes(x=chromvar2_RORC.351, y=adt_CD196) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept = 1, color = 'red', linetype ='dashed')+
  geom_hline(yintercept = 1.5, color = 'red', linetype ='dashed')+
  theme_cowplot()+
  labs(x = 'RORC (ATAC)', y = 'CD196 (ADT)', fill = 'Density')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_1/density_17.pdf', width = 4, height = 4)
p1_new
dev.off()

select.cells_17_bi <- rownames(p1_data[p1_data$chromvar2_RORC.351 > 1 & p1_data$adt_CD196 > 1.5,])

p2 <- FeatureScatter(data, 'chromvar2_TBX21-124', 'adt_CD183')
p2_data <- p2$data

p2_new <- ggplot(p2_data, aes(x=chromvar2_TBX21.124, y=adt_CD183) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept = 0.3, color = 'red', linetype ='dashed')+
  geom_hline(yintercept = 1.6, color = 'red', linetype ='dashed')+
  theme_cowplot()+
  labs(x = 'TBX21 (ATAC)', y = 'CD183 (ADT)', fill = 'Density')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_1/density_1.pdf', width = 4, height = 4)
p2_new
dev.off()

select.cells_1_bi <- rownames(p2_data[p2_data$chromvar2_TBX21.124 > 0.3 & p2_data$adt_CD183 > 1.6,])
select.cells_17_1_bi <- intersect(select.cells_17_bi, select.cells_1_bi)

data$celltype_sub <- 'RORC-CD196-TBX21-CD183-'
data$celltype_sub[select.cells_17_bi] <- 'RORC+CD196+TBX21-CD183-'
data$celltype_sub[select.cells_1_bi] <- 'RORC-CD196-TBX21+CD183+'
data$celltype_sub[select.cells_17_1_bi] <- 'RORC+CD196+TBX21+CD183+'
data$celltype_sub <- factor(data$celltype_sub, levels = c('RORC-CD196-TBX21-CD183-', 'RORC-CD196-TBX21+CD183+', 'RORC+CD196+TBX21-CD183-', 'RORC+CD196+TBX21+CD183+'))
Idents(data) <- 'celltype_sub'
pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_1/scatter_RNA.pdf', width = 4, height = 4)
FeatureScatter(data, 'sct_IL17F', 'sct_IFNG', plot.cor = F) + labs(x= 'IL17F (RNA)', y = 'IFNG (RNA)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_1/scatter_ADT.pdf', width = 4, height = 4)
FeatureScatter(data, 'adt_CD196', 'adt_CD183', plot.cor = F) + labs(x= 'CD196 (ADT)', y = 'CD183 (ADT)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_1/scatter_ATAC.pdf', width = 4, height = 4)
FeatureScatter(data, 'chromvar2_RORC-351', 'chromvar2_TBX21-124', plot.cor = F) + labs(x= 'RORC (ATAC)', y = 'TBX21 (ATAC)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_1/UMAP.pdf', width = 4, height = 4)
DimPlot(data, reduction = 'atac.umap', raster = T) + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

for (i in 1:4){
  p <- DimPlot(data, reduction = 'atac.umap', raster = T, cells.highlight = colnames(data)[data$celltype_sub == levels(data$celltype_sub)[i]]) + NoLegend() + ggtitle(levels(data$celltype_sub)[i]) + theme(plot.title = element_text(hjust = 0.5))
  assign(paste0('p', i), p)
}

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_1/UMAP_split.pdf', width = 8, height = 8)
wrap_plots(p1,p2,p3,p4,ncol = 2)
dev.off()

###########
data <- data_bp[,!is.na(data_bp$monocle_Th17_Th1_2)]
rna <- t(data@assays$SCT@data)
adt <- t(data@assays$ADT@data)
motif <- t(data@assays$chromvar2@data)
meta <- data@meta.data

p1 <- plot_marker_condition('monocle_Th17_Th1_2', 'rna', 'IL17F', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p2 <- plot_marker_condition('monocle_Th17_Th1_2', 'rna', 'IL23R', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p3 <- plot_marker_condition('monocle_Th17_Th1_2', 'rna', 'IFNG', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p4 <- plot_marker_condition('monocle_Th17_Th1_2', 'motif', 'RORC', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p5 <- plot_marker_condition('monocle_Th17_Th1_2', 'motif', 'RORA', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'))
p6 <- plot_marker_condition('monocle_Th17_Th1_2', 'motif', 'TBX21', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p7 <- plot_marker_condition('monocle_Th17_Th1_2', 'adt', 'CD196', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p8 <- plot_marker_condition('monocle_Th17_Th1_2', 'adt', 'CD194', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p9 <- plot_marker_condition('monocle_Th17_Th1_2', 'adt', 'CD183', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_2/pathogenic_Th17.pdf', width = 12, height = 12)
wrap_plots(p1, p2, p3, p7, p8, p9, p4, p5, p6, ncol = 3)
dev.off()

p1 <- FeatureScatter(data, 'chromvar2_RORC-351', 'adt_CD196')
p1_data <- p1$data

p1_new <- ggplot(p1_data, aes(x=chromvar2_RORC.351, y=adt_CD196) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept = 1, color = 'red', linetype ='dashed')+
  geom_hline(yintercept = 1.5, color = 'red', linetype ='dashed')+
  theme_cowplot()+
  labs(x = 'RORC (ATAC)', y = 'CD196 (ADT)', fill = 'Density')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_2/density_17.pdf', width = 4, height = 4)
p1_new
dev.off()

select.cells_17_bi <- rownames(p1_data[p1_data$chromvar2_RORC.351 > 1 & p1_data$adt_CD196 > 1.5,])

p2 <- FeatureScatter(data, 'chromvar2_TBX21-124', 'adt_CD183')
p2_data <- p2$data

p2_new <- ggplot(p2_data, aes(x=chromvar2_TBX21.124, y=adt_CD183) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept = 0.3, color = 'red', linetype ='dashed')+
  geom_hline(yintercept = 1.6, color = 'red', linetype ='dashed')+
  theme_cowplot()+
  labs(x = 'TBX21 (ATAC)', y = 'CD183 (ADT)', fill = 'Density')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_2/density_1.pdf', width = 4, height = 4)
p2_new
dev.off()

select.cells_1_bi <- rownames(p2_data[p2_data$chromvar2_TBX21.124 > 0.3 & p2_data$adt_CD183 > 1.6,])
select.cells_17_1_bi <- intersect(select.cells_17_bi, select.cells_1_bi)

data$celltype_sub <- 'RORC-CD196-TBX21-CD183-'
data$celltype_sub[select.cells_17_bi] <- 'RORC+CD196+TBX21-CD183-'
data$celltype_sub[select.cells_1_bi] <- 'RORC-CD196-TBX21+CD183+'
data$celltype_sub[select.cells_17_1_bi] <- 'RORC+CD196+TBX21+CD183+'
data$celltype_sub <- factor(data$celltype_sub, levels = c('RORC-CD196-TBX21-CD183-', 'RORC-CD196-TBX21+CD183+', 'RORC+CD196+TBX21-CD183-', 'RORC+CD196+TBX21+CD183+'))
Idents(data) <- 'celltype_sub'
pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_2/scatter_RNA.pdf', width = 4, height = 4)
FeatureScatter(data, 'sct_IL17F', 'sct_IFNG', plot.cor = F) + labs(x= 'IL17F (RNA)', y = 'IFNG (RNA)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_2/scatter_ADT.pdf', width = 4, height = 4)
FeatureScatter(data, 'adt_CD196', 'adt_CD183', plot.cor = F) + labs(x= 'CD196 (ADT)', y = 'CD183 (ADT)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_2/scatter_ATAC.pdf', width = 4, height = 4)
FeatureScatter(data, 'chromvar2_RORC-351', 'chromvar2_TBX21-124', plot.cor = F) + labs(x= 'RORC (ATAC)', y = 'TBX21 (ATAC)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_2/UMAP.pdf', width = 4, height = 4)
DimPlot(data, reduction = 'atac.umap', raster = T) + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

for (i in 1:4){
  p <- DimPlot(data, reduction = 'atac.umap', raster = T, cells.highlight = colnames(data)[data$celltype_sub == levels(data$celltype_sub)[i]]) + NoLegend() + ggtitle(levels(data$celltype_sub)[i]) + theme(plot.title = element_text(hjust = 0.5))
  assign(paste0('p', i), p)
}

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_2/UMAP_split.pdf', width = 8, height = 8)
wrap_plots(p1,p2,p3,p4,ncol = 2)
dev.off()

###########
data <- data_bp[,!is.na(data_bp$monocle_Th17_Th1_6)]
rna <- t(data@assays$SCT@data)
adt <- t(data@assays$ADT@data)
motif <- t(data@assays$chromvar2@data)
meta <- data@meta.data

p1 <- plot_marker_condition('monocle_Th17_Th1_6', 'rna', 'IL17F', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p2 <- plot_marker_condition('monocle_Th17_Th1_6', 'rna', 'IL23R', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p3 <- plot_marker_condition('monocle_Th17_Th1_6', 'rna', 'IFNG', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p4 <- plot_marker_condition('monocle_Th17_Th1_6', 'motif', 'RORC', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p5 <- plot_marker_condition('monocle_Th17_Th1_6', 'motif', 'RORA', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'))
p6 <- plot_marker_condition('monocle_Th17_Th1_6', 'motif', 'TBX21', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p7 <- plot_marker_condition('monocle_Th17_Th1_6', 'adt', 'CD196', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p8 <- plot_marker_condition('monocle_Th17_Th1_6', 'adt', 'CD194', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p9 <- plot_marker_condition('monocle_Th17_Th1_6', 'adt', 'CD183', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6/pathogenic_Th17.pdf', width = 12, height = 12)
wrap_plots(p1, p2, p3, p7, p8, p9, p4, p5, p6, ncol = 3)
dev.off()

p1 <- FeatureScatter(data, 'chromvar2_RORC-351', 'adt_CD196')
p1_data <- p1$data

p1_new <- ggplot(p1_data, aes(x=chromvar2_RORC.351, y=adt_CD196) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept = 1.6, color = 'red', linetype ='dashed')+
  geom_hline(yintercept = 1.6, color = 'red', linetype ='dashed')+
  theme_cowplot()+
  labs(x = 'RORC (ATAC)', y = 'CD196 (ADT)', fill = 'Density')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6/density_17.pdf', width = 4, height = 4)
p1_new
dev.off()

select.cells_17_bi <- rownames(p1_data[p1_data$chromvar2_RORC.351 > 1.6 & p1_data$adt_CD196 > 1.6,])

p2 <- FeatureScatter(data, 'chromvar2_TBX21-124', 'adt_CD183')
p2_data <- p2$data

p2_new <- ggplot(p2_data, aes(x=chromvar2_TBX21.124, y=adt_CD183) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept = 0.3, color = 'red', linetype ='dashed')+
  geom_hline(yintercept = 1.5, color = 'red', linetype ='dashed')+
  theme_cowplot()+
  labs(x = 'TBX21 (ATAC)', y = 'CD183 (ADT)', fill = 'Density')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6/density_1.pdf', width = 4, height = 4)
p2_new
dev.off()

select.cells_1_bi <- rownames(p2_data[p2_data$chromvar2_TBX21.124 > 0.3 & p2_data$adt_CD183 > 1.5,])
select.cells_17_1_bi <- intersect(select.cells_17_bi, select.cells_1_bi)

data$celltype_sub <- 'RORC-CD196-TBX21-CD183-'
data$celltype_sub[select.cells_17_bi] <- 'RORC+CD196+TBX21-CD183-'
data$celltype_sub[select.cells_1_bi] <- 'RORC-CD196-TBX21+CD183+'
data$celltype_sub[select.cells_17_1_bi] <- 'RORC+CD196+TBX21+CD183+'
data$celltype_sub <- factor(data$celltype_sub, levels = c('RORC-CD196-TBX21-CD183-', 'RORC-CD196-TBX21+CD183+', 'RORC+CD196+TBX21-CD183-', 'RORC+CD196+TBX21+CD183+'))
Idents(data) <- 'celltype_sub'
pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6/scatter_RNA.pdf', width = 4, height = 4)
FeatureScatter(data, 'sct_IL17F', 'sct_IFNG', plot.cor = F) + labs(x= 'IL17F (RNA)', y = 'IFNG (RNA)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6/scatter_ADT.pdf', width = 4, height = 4)
FeatureScatter(data, 'adt_CD196', 'adt_CD183', plot.cor = F) + labs(x= 'CD196 (ADT)', y = 'CD183 (ADT)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6/scatter_ATAC.pdf', width = 4, height = 4)
FeatureScatter(data, 'chromvar2_RORC-351', 'chromvar2_TBX21-124', plot.cor = F) + labs(x= 'RORC (ATAC)', y = 'TBX21 (ATAC)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6/UMAP.pdf', width = 4, height = 4)
DimPlot(data, reduction = 'atac.umap', raster = T) + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

for (i in 1:4){
  p <- DimPlot(data, reduction = 'atac.umap', raster = T, cells.highlight = colnames(data)[data$celltype_sub == levels(data$celltype_sub)[i]]) + NoLegend() + ggtitle(levels(data$celltype_sub)[i]) + theme(plot.title = element_text(hjust = 0.5))
  assign(paste0('p', i), p)
}

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6/UMAP_split.pdf', width = 8, height = 8)
wrap_plots(p1,p2,p3,p4,ncol = 2)
dev.off()

###########
data <- data_bp[,!is.na(data_bp$monocle_Th17_Th1_6)]
data <- data[,data$monocle_Th17_Th1_6 > 10]
rna <- t(data@assays$SCT@data)
adt <- t(data@assays$ADT@data)
motif <- t(data@assays$chromvar2@data)
meta <- data@meta.data

p1 <- plot_marker_condition('monocle_Th17_Th1_6', 'rna', 'IL17F', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p2 <- plot_marker_condition('monocle_Th17_Th1_6', 'rna', 'IL23R', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p3 <- plot_marker_condition('monocle_Th17_Th1_6', 'rna', 'IFNG', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p4 <- plot_marker_condition('monocle_Th17_Th1_6', 'motif', 'RORC', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p5 <- plot_marker_condition('monocle_Th17_Th1_6', 'motif', 'RORA', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'))
p6 <- plot_marker_condition('monocle_Th17_Th1_6', 'motif', 'TBX21', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p7 <- plot_marker_condition('monocle_Th17_Th1_6', 'adt', 'CD196', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p8 <- plot_marker_condition('monocle_Th17_Th1_6', 'adt', 'CD194', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')
p9 <- plot_marker_condition('monocle_Th17_Th1_6', 'adt', 'CD183', c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2')) + theme(legend.position = 'none')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6_1/pathogenic_Th17.pdf', width = 12, height = 12)
wrap_plots(p1, p2, p3, p7, p8, p9, p4, p5, p6, ncol = 3)
dev.off()

p1 <- FeatureScatter(data, 'chromvar2_RORC-351', 'adt_CD196')
p1_data <- p1$data

p1_new <- ggplot(p1_data, aes(x=chromvar2_RORC.351, y=adt_CD196) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept = 0.5, color = 'red', linetype ='dashed')+
  geom_hline(yintercept = 1.53, color = 'red', linetype ='dashed')+
  theme_cowplot()+
  labs(x = 'RORC (ATAC)', y = 'CD196 (ADT)', fill = 'Density')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6_1/density_17.pdf', width = 4, height = 4)
p1_new
dev.off()

select.cells_17_bi <- rownames(p1_data[p1_data$chromvar2_RORC.351 > 0.5 & p1_data$adt_CD196 > 1.53,])

p2 <- FeatureScatter(data, 'chromvar2_TBX21-124', 'adt_CD183')
p2_data <- p2$data

p2_new <- ggplot(p2_data, aes(x=chromvar2_TBX21.124, y=adt_CD183) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept = -1.3, color = 'red', linetype ='dashed')+
  geom_hline(yintercept = 1, color = 'red', linetype ='dashed')+
  theme_cowplot()+
  labs(x = 'TBX21 (ATAC)', y = 'CD183 (ADT)', fill = 'Density')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6_1/density_1.pdf', width = 4, height = 4)
p2_new
dev.off()

select.cells_1_bi <- rownames(p2_data[p2_data$chromvar2_TBX21.124 > -1.3 & p2_data$adt_CD183 > 1,])
select.cells_17_1_bi <- intersect(select.cells_17_bi, select.cells_1_bi)

data$celltype_sub <- 'RORC-CD196-TBX21-CD183-'
data$celltype_sub[select.cells_17_bi] <- 'RORC+CD196+TBX21-CD183-'
data$celltype_sub[select.cells_1_bi] <- 'RORC-CD196-TBX21+CD183+'
data$celltype_sub[select.cells_17_1_bi] <- 'RORC+CD196+TBX21+CD183+'
data$celltype_sub <- factor(data$celltype_sub, levels = c('RORC-CD196-TBX21-CD183-', 'RORC-CD196-TBX21+CD183+', 'RORC+CD196+TBX21-CD183-', 'RORC+CD196+TBX21+CD183+'))
Idents(data) <- 'celltype_sub'
pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6_1/scatter_RNA.pdf', width = 4, height = 4)
FeatureScatter(data, 'sct_IL17F', 'sct_IFNG', plot.cor = F) + labs(x= 'IL17F (RNA)', y = 'IFNG (RNA)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6_1/scatter_ADT.pdf', width = 4, height = 4)
FeatureScatter(data, 'adt_CD196', 'adt_CD183', plot.cor = F) + labs(x= 'CD196 (ADT)', y = 'CD183 (ADT)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6_1/scatter_ATAC.pdf', width = 4, height = 4)
FeatureScatter(data, 'chromvar2_RORC-351', 'chromvar2_TBX21-124', plot.cor = F) + labs(x= 'RORC (ATAC)', y = 'TBX21 (ATAC)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6_1/UMAP.pdf', width = 4, height = 4)
DimPlot(data, reduction = 'atac.umap', raster = F) + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

for (i in 1:4){
  p <- DimPlot(data, reduction = 'atac.umap', raster = F, cells.highlight = colnames(data)[data$celltype_sub == levels(data$celltype_sub)[i]]) + NoLegend() + ggtitle(levels(data$celltype_sub)[i]) + theme(plot.title = element_text(hjust = 0.5))
  assign(paste0('p', i), p)
}

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_6_1/UMAP_split.pdf', width = 8, height = 8)
wrap_plots(p1,p2,p3,p4,ncol = 2)
dev.off()
