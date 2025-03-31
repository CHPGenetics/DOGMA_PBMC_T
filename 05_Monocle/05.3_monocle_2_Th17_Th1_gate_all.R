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

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

data <- readRDS('../output/tcell_annotated_updated_2_conditions.RDS')
add_pseudotime <- function(name){
  pseudotime <- read.csv(paste0('../output/monocle_2/', name, '.csv'), row.names = 'X')
  data <- AddMetaData(data, pseudotime, paste0('monocle_', name))
  return(data)
}

data <- add_pseudotime('cd4_memory')

data <- data[,!is.na(data$monocle_cd4_memory)]

motif <- readRDS('../output/ArchR/DOGMA_filtered_multiome/motif.RDS')
motif <- assays(motif)[[2]]
colnames(motif) <- gsub('#', '_', colnames(motif))
motif <- motif[,colnames(data)]
rownames(motif) <- gsub('_', '-', rownames(motif))
data[['chromvar2']] <- CreateAssayObject(data = motif)
rm(motif)

saveRDS(data, file = '../output/monocle_2/tcell_annotated_updated_2_conditions_cd4_memory.RDS')

p1 <- FeatureScatter(data, 'chromvar2_RORC-351', 'adt_CD196')
p1_data <- p1$data

p1_new <- ggplot(p1_data, aes(x=chromvar2_RORC.351, y=adt_CD196) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept = 1, color = 'red', linetype ='dashed')+
  geom_hline(yintercept = 1.5, color = 'red', linetype ='dashed')+
  theme_cowplot()+
  labs(x = 'RORC (ATAC)', y = 'CD196 (ADT)', fill = 'Density')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC/density_17.pdf', width = 4, height = 4)
p1_new
dev.off()

select.cells_17_bi <- rownames(p1_data[p1_data$chromvar2_RORC.351 > 1 & p1_data$adt_CD196 > 1.5,])

p2 <- FeatureScatter(data, 'chromvar2_TBX21-124', 'adt_CD183')
p2_data <- p2$data

p2_new <- ggplot(p2_data, aes(x=chromvar2_TBX21.124, y=adt_CD183) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept = 0.4, color = 'red', linetype ='dashed')+
  geom_hline(yintercept = 1.6, color = 'red', linetype ='dashed')+
  theme_cowplot()+
  labs(x = 'TBX21 (ATAC)', y = 'CD183 (ADT)', fill = 'Density')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC/density_1.pdf', width = 4, height = 4)
p2_new
dev.off()

select.cells_1_bi <- rownames(p2_data[p2_data$chromvar2_TBX21.124 > 0.4 & p2_data$adt_CD183 > 1.6,])
select.cells_17_1_bi <- intersect(select.cells_17_bi, select.cells_1_bi)

data$celltype_sub <- 'RORC-CD196-TBX21-CD183-'
data$celltype_sub[select.cells_17_bi] <- 'RORC+CD196+TBX21-CD183-'
data$celltype_sub[select.cells_1_bi] <- 'RORC-CD196-TBX21+CD183+'
data$celltype_sub[select.cells_17_1_bi] <- 'RORC+CD196+TBX21+CD183+'
data$celltype_sub <- factor(data$celltype_sub, levels = c('RORC-CD196-TBX21-CD183-', 'RORC-CD196-TBX21+CD183+', 'RORC+CD196+TBX21-CD183-', 'RORC+CD196+TBX21+CD183+'))
Idents(data) <- 'celltype_sub'
pdf('../plots/monocle_2/cd4_memory/ADT_ATAC/scatter_RNA.pdf', width = 4, height = 4)
FeatureScatter(data, 'sct_IL17F', 'sct_IFNG', plot.cor = F) + labs(x= 'IL17F (RNA)', y = 'IFNG (RNA)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC/scatter_ADT.pdf', width = 4, height = 4)
FeatureScatter(data, 'adt_CD196', 'adt_CD183', plot.cor = F) + labs(x= 'CD196 (ADT)', y = 'CD183 (ADT)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC/scatter_ATAC.pdf', width = 4, height = 4)
FeatureScatter(data, 'chromvar2_RORC-351', 'chromvar2_TBX21-124', plot.cor = F) + labs(x= 'RORC (ATAC)', y = 'TBX21 (ATAC)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC/UMAP.pdf', width = 4, height = 4)
DimPlot(data, reduction = 'atac.umap', raster = T) + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

for (i in 1:4){
  p <- DimPlot(data, reduction = 'atac.umap', raster = T, cells.highlight = colnames(data)[data$celltype_sub == levels(data$celltype_sub)[i]]) + NoLegend() + ggtitle(levels(data$celltype_sub)[i]) + theme(plot.title = element_text(hjust = 0.5))
  assign(paste0('p', i), p)
}

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC/UMAP_split.pdf', width = 8, height = 8)
wrap_plots(p1,p2,p3,p4,ncol = 2)
dev.off()

#############
p1 <- FeatureScatter(data, 'chromvar2_RORC-351', 'adt_CD196')
p1_data <- p1$data

p1_new <- ggplot(p1_data, aes(x=chromvar2_RORC.351, y=adt_CD196) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept = 0.5, color = 'red', linetype ='dashed')+
  geom_hline(yintercept = 1.53, color = 'red', linetype ='dashed')+
  theme_cowplot()+
  labs(x = 'RORC (ATAC)', y = 'CD196 (ADT)', fill = 'Density')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re/density_17.pdf', width = 4, height = 4)
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

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re/density_1.pdf', width = 4, height = 4)
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
pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re/scatter_RNA.pdf', width = 4, height = 4)
FeatureScatter(data, 'sct_IL17F', 'sct_IFNG', plot.cor = F) + labs(x= 'IL17F (RNA)', y = 'IFNG (RNA)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re/scatter_ADT.pdf', width = 4, height = 4)
FeatureScatter(data, 'adt_CD196', 'adt_CD183', plot.cor = F) + labs(x= 'CD196 (ADT)', y = 'CD183 (ADT)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re/scatter_ATAC.pdf', width = 4, height = 4)
FeatureScatter(data, 'chromvar2_RORC-351', 'chromvar2_TBX21-124', plot.cor = F) + labs(x= 'RORC (ATAC)', y = 'TBX21 (ATAC)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re/UMAP.pdf', width = 4, height = 4)
DimPlot(data, reduction = 'atac.umap', raster = T) + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

for (i in 1:4){
  p <- DimPlot(data, reduction = 'atac.umap', raster = T, cells.highlight = colnames(data)[data$celltype_sub == levels(data$celltype_sub)[i]]) + NoLegend() + ggtitle(levels(data$celltype_sub)[i]) + theme(plot.title = element_text(hjust = 0.5))
  assign(paste0('p', i), p)
}

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re/UMAP_split.pdf', width = 8, height = 8)
wrap_plots(p1,p2,p3,p4,ncol = 2)
dev.off()

#############
p1 <- FeatureScatter(data, 'chromvar2_RORC-351', 'adt_CD196')
p1_data <- p1$data

p1_new <- ggplot(p1_data, aes(x=chromvar2_RORC.351, y=adt_CD196) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept = 1.1, color = 'red', linetype ='dashed')+
  geom_hline(yintercept = 1.3, color = 'red', linetype ='dashed')+
  theme_cowplot()+
  labs(x = 'RORC (ATAC)', y = 'CD196 (ADT)', fill = 'Density')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re_re/density_17.pdf', width = 4, height = 4)
p1_new
dev.off()

select.cells_17_bi <- rownames(p1_data[p1_data$chromvar2_RORC.351 > 1.1 & p1_data$adt_CD196 > 1.3,])

p2 <- FeatureScatter(data, 'chromvar2_TBX21-124', 'adt_CD183')
p2_data <- p2$data

p2_new <- ggplot(p2_data, aes(x=chromvar2_TBX21.124, y=adt_CD183) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  geom_vline(xintercept = -0.6, color = 'red', linetype ='dashed')+
  geom_hline(yintercept = 1.7, color = 'red', linetype ='dashed')+
  theme_cowplot()+
  labs(x = 'TBX21 (ATAC)', y = 'CD183 (ADT)', fill = 'Density')

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re_re/density_1.pdf', width = 4, height = 4)
p2_new
dev.off()

select.cells_1_bi <- rownames(p2_data[p2_data$chromvar2_TBX21.124 > -0.6 & p2_data$adt_CD183 > 1.7,])
select.cells_17_1_bi <- intersect(select.cells_17_bi, select.cells_1_bi)

data$celltype_sub <- 'RORC-CD196-TBX21-CD183-'
data$celltype_sub[select.cells_17_bi] <- 'RORC+CD196+TBX21-CD183-'
data$celltype_sub[select.cells_1_bi] <- 'RORC-CD196-TBX21+CD183+'
data$celltype_sub[select.cells_17_1_bi] <- 'RORC+CD196+TBX21+CD183+'
data$celltype_sub <- factor(data$celltype_sub, levels = c('RORC-CD196-TBX21-CD183-', 'RORC-CD196-TBX21+CD183+', 'RORC+CD196+TBX21-CD183-', 'RORC+CD196+TBX21+CD183+'))
Idents(data) <- 'celltype_sub'
pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re_re/scatter_RNA.pdf', width = 4, height = 4)
FeatureScatter(data, 'sct_IL17F', 'sct_IFNG', plot.cor = F) + labs(x= 'IL17F (RNA)', y = 'IFNG (RNA)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re_re/scatter_ADT.pdf', width = 4, height = 4)
FeatureScatter(data, 'adt_CD196', 'adt_CD183', plot.cor = F) + labs(x= 'CD196 (ADT)', y = 'CD183 (ADT)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re_re/scatter_ATAC.pdf', width = 4, height = 4)
FeatureScatter(data, 'chromvar2_RORC-351', 'chromvar2_TBX21-124', plot.cor = F) + labs(x= 'RORC (ATAC)', y = 'TBX21 (ATAC)', color = '') + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re_re/UMAP.pdf', width = 4, height = 4)
DimPlot(data, reduction = 'atac.umap', raster = T) + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1))
dev.off()

for (i in 1:4){
  p <- DimPlot(data, reduction = 'atac.umap', raster = T, cells.highlight = colnames(data)[data$celltype_sub == levels(data$celltype_sub)[i]]) + NoLegend() + ggtitle(levels(data$celltype_sub)[i]) + theme(plot.title = element_text(hjust = 0.5))
  assign(paste0('p', i), p)
}

pdf('../plots/monocle_2/cd4_memory/ADT_ATAC_re_re/UMAP_split.pdf', width = 8, height = 8)
wrap_plots(p1,p2,p3,p4,ncol = 2)
dev.off()
