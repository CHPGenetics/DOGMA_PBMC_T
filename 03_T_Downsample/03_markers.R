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

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')

plan('multicore')
options(future.globals.maxSize= 60*1024^3)

data <- readRDS('../output/tcell_annotated_updated.RDS')

for (i in 1:20){
  celltype <- levels(data$celltype_updated)[i]
  assign(paste0('p',i),
         DimPlot(data, reduction = 'wnn2.umap', cells.highlight = colnames(data)[data$celltype_updated == celltype]) + NoLegend() + ggtitle(celltype) + theme(plot.title = element_text(hjust = 0.5)))
}

pdf('../plots/annotation_all_updated/celltype.pdf', width = 20, height = 16)
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20, ncol = 5)
dev.off()

pdf('../plots/annotation_all_updated/RNA_markers.pdf', width = 11, height = 5)
DotPlot(data, assay = 'SCT', c('CD4', 'CCR7',
                               'IL2RA','SLC3A2', 'SLC7A5',
                               'IL7R',
                               'CD8A',
                               'IFNG', 'LTA', 'TNF',
                               'IL17F', 'CCL20', 'RORA',
                               'TNFSF8', 'CD40LG',
                               'GZMB', 'KLRG1', 'LAG3', 'CCL3', 'CCL4', 'CCL5',
                               'TRAV1-2', 'KLRB1', 'SLC4A10', 'DPP4',
                               'TRDV2', 'TRGV9', 'TRDC',
                               'IKZF2', 'FOXP3', 'CTLA4'), col.min = 0, col.max = 1)+
  RotatedAxis()+labs(x = 'RNA', y = '')+theme(axis.text.x = element_text(face = 'italic'))
dev.off()

pdf('../plots/annotation_all_updated/ADT_markers.pdf', width = 10, height = 5)
DotPlot(data, assay = 'ADT', c('CD4', 'CD45RA', 'CD62L',
                               'CD44', 'CD25', 'CD69', 'CD71',
                               'CD127',
                               'CD8',
                               'CD45RO',
                               'CD183',
                               'CD194', 'CD196',
                               'CD185', 'CD200', 'CD279',
                               'KLRG1', 'CD107a',
                               'TCR-Va7.2', 'CD161', 'CD26',
                               'TCR-Vd2',
                               'CD152', 'CD39'), col.min = 0, col.max = 1, cols = c("lightgrey", "darkgreen"))+
  RotatedAxis()+labs(x = 'ADT', y = '')
dev.off()

pdf('../plots/annotation_all_updated/ATAC_markers.pdf', width = 6, height = 5)
DotPlot(data, assay = 'chromvar', c('MA0690.1',
                                    'MA1151.1', 'MA0071.1', 'MA0072.1',
                                    'MA0850.1'), col.min = 0, col.max = 1, cols = c("lightgrey", "darkred"))+
  RotatedAxis()+
  scale_x_discrete(labels= c('TBX21 (MA0690.1)', 'RORC (MA1151.1)', 'RORA (MA0071.1)', 'RORA (MA0072.1)', 'FOXP3 (MA0850.1)'))+labs(x = 'ATAC', y = '')
dev.off()

library(openxlsx)
Idents(data) <- 'celltype_updated'

rna.markers <- FindAllMarkers(data, assay = 'SCT', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.xlsx(rna.markers, '../output/annotation_all_updated/rna_markers_annotated.xlsx', colNames = T, rowNames = F, overwrite = T)
# rna.markers <- read.xlsx('../output/annotation_all_updated/rna_markers_annotated.xlsx')
rna.top5 <- rna.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

pdf('../plots/annotation_all_updated/heatmap_RNA_annotated_downsample.pdf', width = 5, height = 13)
print(DoHeatmap(subset(data, downsample = 700), assay = 'SCT', features = rna.top5$gene, angle = 90, size = 4) + theme(legend.position = 'bottom') + guides(color=FALSE) + theme(axis.text.y = element_text(face = 'italic')))
dev.off()

adt.markers <- FindAllMarkers(data, assay = 'ADT', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.xlsx(adt.markers, '../output/annotation_all_updated/adt_markers_annotated.xlsx', colNames = T, rowNames = F, overwrite = T)
# adt.markers <- read.xlsx('../output/annotation_all_updated/adt_markers_annotated.xlsx')
adt.top5 <- adt.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

pdf('../plots/annotation_all_updated/heatmap_ADT_annotated_downsample.pdf', width = 5, height = 13)
print(DoHeatmap(subset(data, downsample = 700), assay = 'ADT', features = adt.top5$gene, angle = 90, size = 4) + theme(legend.position = 'bottom') + guides(color=FALSE))
dev.off()

motif.markers <- FindAllMarkers(data, assay = 'chromvar', only.pos = TRUE,
                                mean.fxn = rowMeans,
                                fc.name = "avg_diff")
motif.markers$motif <- ConvertMotifID(data, id = motif.markers$gene, assay = 'peaks')
motif.markers$name <- paste0(motif.markers$motif, ' (', motif.markers$gene, ')')
write.xlsx(motif.markers, '../output/annotation_all_updated/motif_markers_annotated.xlsx', colNames = T, rowNames = F, overwrite = T)
# motif.markers <- read.xlsx('../output/annotation_all_updated/motif_markers_annotated.xlsx')
motif.top5 <- motif.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_diff)

motif.markers_uni <- motif.markers[!duplicated(motif.markers$name),]
motif <- data@assays$chromvar@data[motif.markers_uni$gene,]
table(rownames(motif) == motif.markers_uni$gene)
rownames(motif) <- motif.markers_uni$name
data[['chromvar2']] <- CreateAssayObject(data = motif)

pdf('../plots/annotation_all_updated/heatmap_motif_annotated_new_downsample.pdf', width = 5, height = 13)
print(DoHeatmap(subset(data, downsample = 700), assay = 'chromvar2', slot = 'data', features = motif.top5$name, angle = 90, size = 4) + theme(legend.position = 'bottom') + guides(color=FALSE))
dev.off()

p1 <- DimPlot(data, reduction = 'wnn2.umap', label = T, repel = T, label.size = 4)+NoLegend()
pdf('../plots/annotation_all_updated/annotated_updated.pdf', width = 8, height = 8)
print(p1)
dev.off()

p1 <- DimPlot(data, reduction = 'wnn2.umap', label = T, repel = T, label.size = 4) + guides(col = guide_legend(ncol = 1, override.aes = list(size=4)))
pdf('../plots/annotation_all_updated/annotated_updated_legend.pdf', width = 10, height = 8)
print(p1)
dev.off()

table <- read.xlsx('../../Single_Cell_Multimodal_Omics_Experiments_202200721.xlsx')
table <- table[table$Experiment_Name %in% data$orig.ident,]
data$sample <- factor(data$sample, levels = unique(table$Human_Subject_SB_Identifier))

Idents(data) <- 'sample'
p1 <- DimPlot(data, reduction = 'wnn2.umap') + guides(col = guide_legend(ncol = 1, override.aes = list(size=4)))
Idents(data) <- 'condition'
p2 <- DimPlot(data, reduction = 'wnn2.umap')
pdf('../plots/annotation_all_updated/annotated_updated_sample_condition.pdf', width = 20, height = 8)
print(p1 | p2)
dev.off()

Idents(data) <- 'sample'
p1 <- VlnPlot(data, 'nFeature_RNA', pt.size = 0, log = T) + NoLegend() + ggtitle('') + xlab('') + ylab('# of genes') + coord_flip()
p2 <- VlnPlot(data, 'nCount_RNA', pt.size = 0, log = T) + NoLegend() + ggtitle('') + xlab('') + ylab('# of UMIs') + coord_flip()
p3 <- VlnPlot(data, 'percent.mt', pt.size = 0) + NoLegend() + ggtitle('') + xlab('') + ylab('%UMIs from mtRNA') + coord_flip()
p4 <- VlnPlot(data, 'sum.ctrl', pt.size = 0) + NoLegend() + ggtitle('') + xlab('') + ylab('# of tags from isotype controls') + coord_flip()
p5 <- VlnPlot(data, 'atac_fragments', pt.size = 0, log = T) + NoLegend() + ggtitle('') + xlab('') + ylab('# of fragments') + coord_flip()
p6 <- VlnPlot(data, 'TSS.enrichment', pt.size = 0, log = T) + NoLegend() + ggtitle('') + xlab('') + ylab('TSS enrichment score') + coord_flip()
pdf('../plots/annotation_all_updated/violin_QC.pdf', width = 18, height = 5)
print(wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 6))
dev.off()

Idents(data) <- 'celltype_updated'
p1 <- VlnPlot(data, 'Activated1', pt.size = 0) + ggtitle('Activated') + NoLegend() + xlab('') + coord_flip()
p2 <- VlnPlot(data, 'Resting1', pt.size = 0) + ggtitle('Resting') + NoLegend() + xlab('') + coord_flip()
pdf('../plots/annotation_all_updated/vln_activated_resting.pdf', width = 8, height = 8)
print(p1 / p2)
dev.off()

p1 <- FeaturePlot(data, 'Activated1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle('Activated')
p2 <- FeaturePlot(data, 'Resting1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle('Resting') 
pdf(paste0('../plots/annotation_all_updated/activated_resting.pdf'), width = 8, height = 4)
print(p1 | p2)
dev.off()

check_marker <- function(type, marker, name = str_split_fixed(marker, '-', 2)[,1]){
  if (type == 'rna'){
    p <- FeaturePlot(data, paste0('sct_', marker), reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle(paste0(name, ' (RNA)'))
  }
  if (type == 'adt'){
    p <- FeaturePlot(data, paste0('adt_', marker), reduction = 'wnn2.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle(paste0(name, ' (ADT)'))
  }
  if (type == 'motif'){
    p <- FeaturePlot(data, paste0('chromvar_', marker), reduction = 'wnn2.umap', cols = c("lightgrey", "darkred"), min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle(paste0(name, ' (ATAC)'))
  }
  return(p)
}

pdf(paste0('../plots/annotation_all_updated/cd25_cd127.pdf'), width = 8, height = 4)
print(wrap_plots(check_marker('adt', 'CD25'), check_marker('adt', 'CD127'), ncol = 2))
dev.off()

data_plot <- table(data$celltype_updated, data$sample)
data_plot <- reshape2::melt(prop.table(data_plot, 2))
colnames(data_plot) <- c('Clusters', 'Donors', 'Proportions')
data_plot$Proportions <- data_plot$Proportions*100

p1 <- ggplot(data_plot, aes(x = Proportions, y = Donors, fill = Clusters))+
  geom_col(position = 'fill', col = 'white')+
  theme_bw()
pdf('../plots/annotation_all_updated/donor_proportion_annotated_updated.pdf', width = 10, height = 5)
print(p1)
dev.off()

data_plot <- table(data$condition, data$sample)
data_plot <- reshape2::melt(prop.table(data_plot, 2))
colnames(data_plot) <- c('Conditions', 'Donors', 'Proportions')
data_plot$Proportions <- data_plot$Proportions*100

p1 <- ggplot(data_plot, aes(x = Proportions, y = Donors, fill = Conditions))+
  geom_col(position = 'fill', col = 'white')+
  theme_bw()
pdf('../plots/annotation_all_updated/donor_proportion_condition_updated.pdf', width = 10, height = 5)
print(p1)
dev.off()

Idents(data) <- 'celltype_updated'
p1 <- DimPlot(data, reduction = 'wnn.umap', label = T, repel = T, label.size = 4) + NoLegend() + ggtitle('WNN') + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(data, reduction = 'rna.umap', label = T, repel = T, label.size = 4) + NoLegend() + ggtitle('RNA') + theme(plot.title = element_text(hjust = 0.5))
p3 <- DimPlot(data, reduction = 'adt.umap', label = T, repel = T, label.size = 4) + NoLegend() + ggtitle('ADT') + theme(plot.title = element_text(hjust = 0.5))
p4 <- DimPlot(data, reduction = 'atac.umap', label = T, repel = T, label.size = 4) + ggtitle('ATAC') + theme(plot.title = element_text(hjust = 0.5))
pdf('../plots/annotation_all_updated/annotated_updated_all.pdf', width = 34, height = 8)
print(wrap_plots(p1, p2, p3, p4, ncol = 4))
dev.off()