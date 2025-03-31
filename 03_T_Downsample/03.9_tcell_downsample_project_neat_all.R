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

data <- readRDS('../output/tcell_annotated_updated_downsample.RDS')

data <- RunUMAP(data, nn.name = "weighted.nn2", reduction.name = "wnn2.umap", reduction.key = "wnn2UMAP_", return.model = T)

neat <- readRDS('~/RWorkSpace/GEO_data/NEAT-seq/Analysis/data_peak_calling_updated.RDS')

anchors <- FindTransferAnchors(
    reference = data,
    query = neat,
    reference.reduction = "lsi",
    reduction = 'lsiproject',
    dims = 2:30,
    reference.assay = 'peaks',
    query.assay = 'peaks'
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

saveRDS(data, file = '../output/downsample/projected_neat_all.RDS')

library(scales)
neat$predicted.id=factor(neat$predicted.id, levels = 0:10)
DimPlot(neat, reduction = "ref.umap", group.by = "predicted.id", cols = hue_pal()(11)[c(1,2,4,6,7)])
DimPlot(data, reduction = "wnn2.umap")

xlim=range(data@reductions$wnn2.umap@cell.embeddings[,1])
ylim=range(data@reductions$wnn2.umap@cell.embeddings[,2])
pdf('../plots/downsample/project_neat/UMAP_all.pdf', width = 17, height = 5)
Idents(data)='harmony_wnn2_clusters'
Idents(neat)="predicted.id"
p2=DimPlot(data, reduction = 'wnn2.umap', raster = T) + xlim(xlim) + ylim(ylim) + NoLegend() + ggtitle('Cells') + theme(plot.title = element_text(hjust = 0.5))
p3=DimPlot(neat, reduction = 'ref.umap', raster = T, cols = hue_pal()(11)[c(1,2,4,6,7)]) + xlim(xlim) + ylim(ylim) + NoLegend() + ggtitle('Nuclei') + theme(plot.title = element_text(hjust = 0.5))
Idents(data)='celltype_updated'
p1=DimPlot(data, reduction = 'wnn2.umap', raster = T)
p1|p2|p3
dev.off()

#######
neat <- MapQuery(
  anchorset = anchors,
  query = neat,
  reference = data,
  refdata = data$celltype_updated,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = "wnn2.umap"
)

pdf('../plots/downsample/project_neat/UMAP_celltype_all.pdf', width = 15, height = 5)
neat$predicted.id=factor(neat$predicted.id, levels = levels(data$celltype_updated))
Idents(neat)="predicted.id"
p3=DimPlot(neat, reduction = 'ref.umap', raster = T, cols = hue_pal()(20)[c(1:8,11:13,16:18)]) + xlim(xlim) + ylim(ylim)
Idents(data)='celltype_updated'
p1=DimPlot(data, reduction = 'wnn2.umap', raster = T) + xlim(xlim) + ylim(ylim)
p1|p3
dev.off()
