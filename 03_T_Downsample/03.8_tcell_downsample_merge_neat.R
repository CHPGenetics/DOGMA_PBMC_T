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

data$tissue <- 'Cells'

neat <- readRDS('~/RWorkSpace/GEO_data/NEAT-seq/Analysis/data_peak_calling_updated.RDS')
neat$tissue <- 'Nuclei'
neat$sample <- 'NEAT'

data <- merge(data, neat)
rm(neat)

# DefaultAssay(data) <- 'RNA'
# data <- SCTransform(data, method = "glmGamPoi", vars.to.regress = "percent.mt") %>% RunPCA(reduction.name = "pca") %>%
#   RunHarmony(group.by.vars = c("sample"), max.iter.harmony = 20, reduction = 'pca', assay.use = 'SCT',project.dim = FALSE,  reduction.save = "harmony_RNA")
# data <- RunUMAP(data, reduction = "harmony_RNA", dims = 1:30, assay = 'SCT', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
# data <- FindNeighbors(data, dims = 1:30, reduction = 'harmony_RNA')
# data <- FindClusters(data, resolution = 0.2, graph.name = 'SCT_snn', algorithm = 3)
# data$merged_rna_clusters <- Idents(data)
# 
# DefaultAssay(data) <- 'ADT'
# # # we will use all ADT features for dimensional reduction
# # # we set a dimensional reduction name to avoid overwriting the
# VariableFeatures(data) <- rownames(data[["ADT"]])
# data <- NormalizeData(data, normalization.method = 'CLR', margin = 2) %>%
#   ScaleData() %>% RunPCA(reduction.name = 'apca') %>%
#   RunHarmony(group.by.vars = c("sample", "cryopreservation"), max.iter.harmony = 20, reduction = 'apca', assay.use = 'ADT',project.dim = FALSE,  reduction.save = "harmony_ADT")
# data <- RunUMAP(data, reduction = "harmony_ADT", dims = 1:30, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
# data <- FindNeighbors(data, dims = 1:30, reduction = 'harmony_ADT')
# data <- FindClusters(data, resolution = 0.2, graph.name = 'ADT_snn', algorithm = 3)
# data$merged_adt_clusters <- Idents(data)

DefaultAssay(data) <- "peaks"
data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = 'q75')
data <- RunSVD(data)
data <- RunHarmony(data, group.by.vars = c("tissue"), max.iter.harmony = 20, reduction = 'lsi', assay.use = 'peaks',project.dim = FALSE,  reduction.save = "harmony_peaks")
data <- RunUMAP(data, reduction = "harmony_peaks", dims = 2:30, assay = 'peaks', reduction.name = 'atac.umap', reduction.key = 'atacUMAP_')
data <- FindNeighbors(data, dims = 2:30, reduction = 'harmony_peaks')
data <- FindClusters(data, resolution = 0.2, graph.name = 'peaks_snn', algorithm = 3)
data$merged_atac_clusters <- Idents(data)

# # Now run multimodal neighbors and embedding
# data <- FindMultiModalNeighbors(object = data,
#                                 reduction.list = list("harmony_RNA", "harmony_peaks", "harmony_ADT"),
#                                 dims.list = list(1:30, 2:30, 1:30))
# data <- RunUMAP(data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_" )
# data <- FindClusters(data, graph.name = "wsnn", resolution = 0.2, algorithm = 3)
# 
# data$merged_wnn_clusters <- Idents(data)
# 
# data <- FindMultiModalNeighbors(object = data,
#                                 reduction.list = list("harmony_peaks", "harmony_ADT"),
#                                 dims.list = list(2:30, 1:30),
#                                 knn.graph.name = "wknn2",
#                                 snn.graph.name = "wsnn2",
#                                 weighted.nn.name = "weighted.nn2")
# data <- RunUMAP(data, nn.name = "weighted.nn2", reduction.name = "wnn2.umap", reduction.key = "wnn2UMAP_" )
# data <- FindClusters(data, graph.name = "wsnn2", algorithm = 3, resolution = 0.2)
# data$merged_wnn2_clusters <- Idents(data)

saveRDS(data, file = '../output/downsample/data_downsample_merge_neat.RDS')

data$celltype_updated = factor(data$celltype_updated, levels =  c('CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
                                                  'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1',
                                                  'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17',
                                                  'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh',
                                                  'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other'))

pdf('../plots/downsample/merge_neat/UMAP.pdf', width = 17, height = 5)
Idents(data)='merged_atac_clusters'
p2=DimPlot(data, reduction = 'atac.umap', cells = colnames(data)[data$tissue == 'Cells'], raster = T) + NoLegend() + ggtitle('Cells') + theme(plot.title = element_text(hjust = 0.5))
p3=DimPlot(data, reduction = 'atac.umap', cells = colnames(data)[data$tissue == 'Nuclei'], raster = T) + NoLegend() + ggtitle('Nuclei') + theme(plot.title = element_text(hjust = 0.5))
Idents(data)='celltype_updated'
p1=DimPlot(data, reduction = 'atac.umap', raster = T, cells = colnames(data)[data$tissue == 'Cells'])
p1|p2|p3
dev.off()
  
pdf('../plots/downsample/merge_neat/activated_resting.pdf', width = 5, height = 10)
p3 <- FeaturePlot(data, 'Activated1', cells = colnames(data)[data$tissue == 'Cells'], reduction = 'atac.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + NoLegend() + ggtitle('Activated')
p4 <- FeaturePlot(data, 'Resting1', cells = colnames(data)[data$tissue == 'Cells'], reduction = 'atac.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + NoLegend() + ggtitle('Resting')
p3/p4
dev.off()

check_marker <- function(type, marker, name = str_split_fixed(marker, '-', 2)[,1], celltype = 'Nuclei'){
  if (type == 'rna'){
    p <- FeaturePlot(data, paste0('sct_', marker), reduction = 'atac.umap', min.cutoff = 'q1', max.cutoff = 'q99', cells = colnames(data)[data$tissue == celltype], raster = T) + ggtitle(paste0(name, ' (RNA)'))
  }
  if (type == 'adt'){
    p <- FeaturePlot(data, paste0('adt_', marker), reduction = 'atac.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q1', max.cutoff = 'q99', cells = colnames(data)[data$tissue == celltype], raster = T) + ggtitle(paste0(name, ' (ADT)'))
  }
  if (type == 'motif'){
    p <- FeaturePlot(data, paste0('chromvar_', marker), reduction = 'atac.umap', cols = c("lightgrey", "darkred"), min.cutoff = 'q1', max.cutoff = 'q99', cells = colnames(data)[data$tissue == celltype], raster = T) + ggtitle(paste0(name, ' (ATAC)'))
  }
  return(p)
}

pdf('../plots/downsample/merge_neat/markers_nuclei.pdf', width = 12, height = 20)
wrap_plots(check_marker('rna', 'RORC'), check_marker('rna', 'GATA3'), check_marker('rna', 'FOXP3'), 
           check_marker('rna', 'TBX21'), check_marker('rna', 'IKZF2'), 
           check_marker('adt', 'RORgT-ADT', 'RORgT'), check_marker('adt', 'GATA3-ADT', 'GATA3'), check_marker('adt', 'FOXP3-ADT', 'FOXP3'), 
           check_marker('adt', 'Tbet-ADT', 'T-bet'), check_marker('adt', 'Helios-ADT', 'Helios'), 
           check_marker('motif', 'MA1151.1', 'RORC'), check_marker('motif', 'MA0037.3', 'GATA3'), check_marker('motif', 'MA0850.1', 'FOXP3'), 
           check_marker('motif', 'MA0690.1', 'TBX21'), check_marker('motif', 'MA1508.1', 'IKZF1'),
           ncol = 3, byrow = F)
dev.off()

check_marker <- function(type, marker, name = str_split_fixed(marker, '-', 2)[,1], celltype = 'Cells'){
  if (type == 'rna'){
    p <- FeaturePlot(data, paste0('sct_', marker), reduction = 'atac.umap', min.cutoff = 'q1', max.cutoff = 'q99', cells = colnames(data)[data$tissue == celltype], raster = T) + ggtitle(paste0(name, ' (RNA)'))
  }
  if (type == 'adt'){
    p <- FeaturePlot(data, paste0('adt_', marker), reduction = 'atac.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q1', max.cutoff = 'q99', cells = colnames(data)[data$tissue == celltype], raster = T) + ggtitle(paste0(name, ' (ADT)'))
  }
  if (type == 'motif'){
    p <- FeaturePlot(data, paste0('chromvar_', marker), reduction = 'atac.umap', cols = c("lightgrey", "darkred"), min.cutoff = 'q1', max.cutoff = 'q99', cells = colnames(data)[data$tissue == celltype], raster = T) + ggtitle(paste0(name, ' (ATAC)'))
  }
  return(p)
}

pdf('../plots/downsample/merge_neat/markers_cells.pdf', width = 12, height = 20)
wrap_plots(check_marker('rna', 'RORC'), check_marker('rna', 'GATA3'), check_marker('rna', 'FOXP3'), 
           check_marker('rna', 'TBX21'), check_marker('rna', 'IKZF2'), 
           check_marker('adt', 'CD196'), check_marker('adt', 'CD194'), check_marker('adt', 'CD152'), 
           check_marker('adt', 'CD183'), check_marker('adt', 'CD39'), 
           check_marker('motif', 'MA1151.1', 'RORC'), check_marker('motif', 'MA0037.3', 'GATA3'), check_marker('motif', 'MA0850.1', 'FOXP3'), 
           check_marker('motif', 'MA0690.1', 'TBX21'), check_marker('motif', 'MA1508.1', 'IKZF1'),
           ncol = 3, byrow = F)
dev.off()
