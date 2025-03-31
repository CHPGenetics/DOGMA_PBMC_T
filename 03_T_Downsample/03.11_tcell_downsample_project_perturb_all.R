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

load('~/RWorkSpace/GEO_data/DOGMA-seq/asap_reproducibility/CD4_CRISPR_asapseq/output/ArchR_subset/Seurat.RData')
neat=data
rm(data)

data <- readRDS('../output/tcell_annotated_updated_downsample.RDS')

data <- RunUMAP(data, nn.name = "weighted.nn2", reduction.name = "wnn2.umap", reduction.key = "wnn2UMAP_", return.model = T)

anchors <- FindTransferAnchors(
    reference = data,
    query = neat,
    reference.reduction = "lsi",
    reduction = 'lsiproject',
    dims = 2:30,
    reference.assay = 'peaks',
    query.assay = 'peaks2'
  )

neat <- MapQuery(
  anchorset = anchors,
  query = neat,
  reference = data,
  refdata = data$harmony_wnn2_clusters,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = "wnn2.umap"
)

saveRDS(neat, file = '../output/downsample/projected_perturb_all.RDS')

library(scales)
neat$predicted.id=factor(neat$predicted.id, levels = 0:10)
DimPlot(neat, reduction = "ref.umap", group.by = "predicted.id", cols = hue_pal()(11)[c(1:5,7,8,11)])
DimPlot(data, reduction = "wnn2.umap")

xlim=range(data@reductions$wnn2.umap@cell.embeddings[,1])
ylim=range(data@reductions$wnn2.umap@cell.embeddings[,2])
pdf('../plots/downsample/project_perturb/UMAP_all.pdf', width = 17, height = 5)
Idents(data)='harmony_wnn2_clusters'
Idents(neat)="predicted.id"
p2=DimPlot(data, reduction = 'wnn2.umap', raster = T) + xlim(xlim) + ylim(ylim) + NoLegend() + ggtitle('In-house') + theme(plot.title = element_text(hjust = 0.5))
p3=DimPlot(neat, reduction = 'ref.umap', raster = T, cols = hue_pal()(11)[c(1:5,7,8,11)]) + xlim(xlim) + ylim(ylim) + NoLegend() + ggtitle('CRISPR') + theme(plot.title = element_text(hjust = 0.5))
Idents(data)='celltype_updated'
p1=DimPlot(data, reduction = 'wnn2.umap', raster = T)
p1|p2|p3
dev.off()

pdf('../plots/downsample/project_perturb/UMAP_all_updated.pdf', width = 16, height = 4)
Idents(data)='celltype_updated'
p1=DimPlot(data, reduction = 'wnn2.umap', raster = T) + xlim(xlim) + ylim(ylim) + guides(col = guide_legend(ncol = 2, override.aes = list(size=4))) + ggtitle('In-house') + theme(plot.title = element_text(hjust = 0.5))
neat$genotype=factor(neat$genotype, levels = paste0('Genotype_', 0:2))
Idents(neat)="genotype"
p2=DimPlot(neat, reduction = 'ref.umap', raster = T) + xlim(xlim) + ylim(ylim) + ggtitle('CRISPR') + theme(plot.title = element_text(hjust = 0.5)) + geom_density_2d(aes(x = wnn2UMAP_1, y = wnn2UMAP_2), data = p1$data, bins = 5)
p1|p2
dev.off()

#######
neat$HTO=gsub('sgGuide',  '', neat$HTO)
neat$HTO=factor(neat$HTO, levels = c('sgNTC_1', 'sgNTC_2',
'sgCD4_1', 'sgCD4_2',
'sgNFKB2_1', 'sgNFKB2_2', 
'sgCD3ECD4_1', 'sgCD3E_2',
'sgZAP70_1', 'sgZAP70_2'))

pdf('../plots/downsample/project_perturb/activated_resting.pdf', width = 16, height = 4)
p1 <- FeaturePlot(data, 'Activated1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + xlim(xlim) + ylim(ylim) + NoLegend() + ggtitle('Activated')
p2 <- FeaturePlot(data, 'Resting1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + xlim(xlim) + ylim(ylim) + NoLegend() + ggtitle('Resting')
Idents(neat) = 'HTO'
p3=DimPlot(neat, reduction = 'ref.umap', raster = T, cells = colnames(neat)[neat$HTO_sgRNA %in% c('sgNTC', 'sgCD4', 'sgNFKB2')]) + xlim(xlim) + ylim(ylim) + ggtitle('sgRNA') + theme(plot.title = element_text(hjust = 0.5))
p4=DimPlot(neat, reduction = 'ref.umap', raster = T, cells = colnames(neat)[!neat$HTO_sgRNA %in% c('sgNTC', 'sgCD4', 'sgNFKB2')]) + xlim(xlim) + ylim(ylim) + ggtitle('sgRNA') + theme(plot.title = element_text(hjust = 0.5))
wrap_plots(p1,p2,p3,p4,ncol = 4)
dev.off()

pdf('../plots/downsample/project_perturb/activated_resting_lined.pdf', width = 16, height = 4)
p5=p1 + geom_density_2d(aes(x = wnn2UMAP_1, y = wnn2UMAP_2), data = p1$data[p1$data$Activated1 > median(p1$data$Activated1),], bins = 5)
p6=p2 + geom_density_2d(aes(x = wnn2UMAP_1, y = wnn2UMAP_2), data = p2$data[p2$data$Resting1 > median(p2$data$Resting1),], bins = 5)
p7=p3 + geom_density_2d(aes(x = wnn2UMAP_1, y = wnn2UMAP_2), data = p1$data[p1$data$Activated1 > median(p1$data$Activated1),], bins = 5)
p8=p4 + geom_density_2d(aes(x = wnn2UMAP_1, y = wnn2UMAP_2), data = p2$data[p2$data$Resting1 > median(p2$data$Resting1),], bins = 5)
wrap_plots(p5,p6,p7,p8,ncol = 4)
dev.off()