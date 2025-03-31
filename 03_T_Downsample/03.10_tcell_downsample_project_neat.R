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

data <- data[,data$celltype_updated %in% c('CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
                                                  'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1',
                                                  'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17',
                                                  'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh',
                                                  'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other')]

DefaultAssay(data) <- "RNA"
data <- SCTransform(data, method = "glmGamPoi", vars.to.regress = "percent.mt") %>% RunPCA(reduction.name = "pca") %>%
  RunHarmony(group.by.vars = c("sample"), max.iter.harmony = 20, reduction = 'pca', assay.use = 'SCT',project.dim = FALSE,  reduction.save = "harmony_RNA")
data <- RunUMAP(data, reduction = "harmony_RNA", dims = 1:30, assay = 'SCT', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
data <- FindNeighbors(data, dims = 1:30, reduction = 'harmony_RNA')
data <- FindClusters(data, resolution = 0.2, graph.name = 'SCT_snn', algorithm = 3)
data$harmony_rna_clusters <- Idents(data)

DefaultAssay(data) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the
VariableFeatures(data) <- rownames(data[["ADT"]])
data <- NormalizeData(data, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca') %>%
  RunHarmony(group.by.vars = c("sample", "cryopreservation"), max.iter.harmony = 20, reduction = 'apca', assay.use = 'ADT',project.dim = FALSE,  reduction.save = "harmony_ADT")
data <- RunUMAP(data, reduction = "harmony_ADT", dims = 1:30, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
data <- FindNeighbors(data, dims = 1:30, reduction = 'harmony_ADT')
data <- FindClusters(data, resolution = 0.2, graph.name = 'ADT_snn', algorithm = 3)
data$harmony_adt_clusters <- Idents(data)

DefaultAssay(data) <- "peaks"
data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = 'q0')
data <- RunSVD(data)
data <- RunHarmony(data, group.by.vars = c("sample"), max.iter.harmony = 20, reduction = 'lsi', assay.use = 'peaks',project.dim = FALSE,  reduction.save = "harmony_peaks")
data <- RunUMAP(data, reduction = "harmony_peaks", dims = 2:30, assay = 'peaks', reduction.name = 'atac.umap', reduction.key = 'atacUMAP_')
data <- FindNeighbors(data, dims = 2:30, reduction = 'harmony_peaks')
data <- FindClusters(data, resolution = 0.2, graph.name = 'peaks_snn', algorithm = 3)
data$harmony_atac_clusters <- Idents(data)

# Now run multimodal neighbors and embedding
data <- FindMultiModalNeighbors(object = data,
                                reduction.list = list("harmony_RNA", "harmony_peaks", "harmony_ADT"),
                                dims.list = list(1:30, 2:30, 1:30))
data <- RunUMAP(data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_" )
data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = 0.2)

data$harmony_wnn_clusters <- Idents(data)

data <- FindMultiModalNeighbors(object = data,
                                reduction.list = list("harmony_peaks", "harmony_ADT"),
                                dims.list = list(2:30, 1:30),
                                knn.graph.name = "wknn2",
                                snn.graph.name = "wsnn2",
                                weighted.nn.name = "weighted.nn2")
data <- RunUMAP(data, nn.name = "weighted.nn2", reduction.name = "wnn2.umap", reduction.key = "wnn2UMAP_" )
data <- FindClusters(data, graph.name = "wsnn2", algorithm = 3, resolution = 0.2)
data$harmony_wnn2_clusters <- Idents(data)

dataset <- 'downsample/project_neat'

p1 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = 'harmony_wnn_clusters')
p2 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "celltype_updated")
pdf(paste0('../plots/', dataset, '/predicted_WNN.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = 'harmony_rna_clusters')
p2 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = "celltype_updated")
pdf(paste0('../plots/', dataset, '/predicted_RNA.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = 'harmony_adt_clusters')
p2 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = "celltype_updated")
pdf(paste0('../plots/', dataset, '/predicted_ADT.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = 'harmony_atac_clusters')
p2 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = "celltype_updated")
pdf(paste0('../plots/', dataset, '/predicted_ATAC.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = 'harmony_wnn2_clusters')
p2 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = "celltype_updated")
pdf(paste0('../plots/', dataset, '/predicted_WNN2.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

saveRDS(data, file = '../output/downsample/data_downsample_project_neat.RDS')

data = readRDS('../output/downsample/data_downsample_project_neat.RDS')

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

saveRDS(data, file = '../output/downsample/projected_neat.RDS')

library(scales)
neat$predicted.id=factor(neat$predicted.id, levels = 0:6)
DimPlot(neat, reduction = "ref.umap", group.by = "predicted.id", cols = hue_pal()(7)[c(1,2,5,6)])
DimPlot(data, reduction = "wnn2.umap")

xlim=range(data@reductions$wnn2.umap@cell.embeddings[,1])
ylim=range(data@reductions$wnn2.umap@cell.embeddings[,2])
pdf('../plots/downsample/project_neat/UMAP.pdf', width = 17, height = 5)
Idents(data)='harmony_wnn2_clusters'
Idents(neat)="predicted.id"
p2=DimPlot(data, reduction = 'wnn2.umap', raster = T) + xlim(xlim) + ylim(ylim) + NoLegend() + ggtitle('Cells') + theme(plot.title = element_text(hjust = 0.5))
p3=DimPlot(neat, reduction = 'ref.umap', raster = T, cols = hue_pal()(7)[c(1,2,5,6)]) + xlim(xlim) + ylim(ylim) + NoLegend() + ggtitle('Nuclei') + theme(plot.title = element_text(hjust = 0.5))
Idents(data)='celltype_updated'
p1=DimPlot(data, reduction = 'wnn2.umap', raster = T)
p1|p2|p3
dev.off()
  
pdf('../plots/downsample/project_neat/activated_resting.pdf', width = 5, height = 10)
p3 <- FeaturePlot(data, 'Activated1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + NoLegend() + ggtitle('Activated')
p4 <- FeaturePlot(data, 'Resting1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + NoLegend() + ggtitle('Resting')
p3/p4
dev.off()

check_marker <- function(type, marker, name = str_split_fixed(marker, '-', 2)[,1]){
  if (type == 'rna'){
    p <- FeaturePlot(neat, paste0('sct_', marker), reduction = 'ref.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + ggtitle(paste0(name, ' (RNA)'))
  }
  if (type == 'adt'){
    p <- FeaturePlot(neat, paste0('adt_', marker), reduction = 'ref.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + ggtitle(paste0(name, ' (ADT)'))
  }
  if (type == 'motif'){
    p <- FeaturePlot(neat, paste0('chromvar_', marker), reduction = 'ref.umap', cols = c("lightgrey", "darkred"), min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + ggtitle(paste0(name, ' (ATAC)'))
  }
  return(p)
}

pdf('../plots/downsample/project_neat/markers_nuclei.pdf', width = 12, height = 20)
wrap_plots(check_marker('rna', 'RORC'), check_marker('rna', 'GATA3'), check_marker('rna', 'FOXP3'), 
           check_marker('rna', 'TBX21'), check_marker('rna', 'IKZF2'), 
           check_marker('adt', 'RORgT-ADT', 'RORgT'), check_marker('adt', 'GATA3-ADT', 'GATA3'), check_marker('adt', 'FOXP3-ADT', 'FOXP3'), 
           check_marker('adt', 'Tbet-ADT', 'T-bet'), check_marker('adt', 'Helios-ADT', 'Helios'), 
           check_marker('motif', 'MA1151.1', 'RORC'), check_marker('motif', 'MA0037.3', 'GATA3'), check_marker('motif', 'MA0850.1', 'FOXP3'), 
           check_marker('motif', 'MA0690.1', 'TBX21'), check_marker('motif', 'MA1508.1', 'IKZF1'),
           ncol = 3, byrow = F)
dev.off()

check_marker <- function(type, marker, name = str_split_fixed(marker, '-', 2)[,1]){
  if (type == 'rna'){
    p <- FeaturePlot(data, paste0('sct_', marker), reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + ggtitle(paste0(name, ' (RNA)'))
  }
  if (type == 'adt'){
    p <- FeaturePlot(data, paste0('adt_', marker), reduction = 'wnn2.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + ggtitle(paste0(name, ' (ADT)'))
  }
  if (type == 'motif'){
    p <- FeaturePlot(data, paste0('chromvar_', marker), reduction = 'wnn2.umap', cols = c("lightgrey", "darkred"), min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + ggtitle(paste0(name, ' (ATAC)'))
  }
  return(p)
}

pdf('../plots/downsample/project_neat/markers_cells.pdf', width = 12, height = 20)
wrap_plots(check_marker('rna', 'RORC'), check_marker('rna', 'GATA3'), check_marker('rna', 'FOXP3'), 
           check_marker('rna', 'TBX21'), check_marker('rna', 'IKZF2'), 
           check_marker('adt', 'CD196'), check_marker('adt', 'CD194'), check_marker('adt', 'CD152'), 
           check_marker('adt', 'CD183'), check_marker('adt', 'CD39'), 
           check_marker('motif', 'MA1151.1', 'RORC'), check_marker('motif', 'MA0037.3', 'GATA3'), check_marker('motif', 'MA0850.1', 'FOXP3'), 
           check_marker('motif', 'MA0690.1', 'TBX21'), check_marker('motif', 'MA1508.1', 'IKZF1'),
           ncol = 3, byrow = F)
dev.off()

########
neat <- MapQuery(
  anchorset = anchors,
  query = neat,
  reference = data,
  refdata = data$celltype_updated,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = "wnn2.umap"
)

pdf('../plots/downsample/project_neat/UMAP_celltype.pdf', width = 15, height = 5)
neat$predicted.id=factor(neat$predicted.id, levels = c('CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
                                                  'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1',
                                                  'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17',
                                                  'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh',
                                                  'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other'))
Idents(neat)="predicted.id"
p3=DimPlot(neat, reduction = 'ref.umap', raster = T, cols = hue_pal()(10)[c(1,2,3,5,9,10)]) + xlim(xlim) + ylim(ylim)
Idents(data)='celltype_updated'
p1=DimPlot(data, reduction = 'wnn2.umap', raster = T) + xlim(xlim) + ylim(ylim)
p1|p3
dev.off()
