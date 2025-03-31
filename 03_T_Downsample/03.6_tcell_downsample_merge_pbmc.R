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

read.marker <- function(path, sheet){
  marker <- read.xlsx(xlsxFile = path, sheet = sheet)
  colnames(marker) <- marker[3,]
  marker <- marker[-c(1:3),]
  marker <- marker %>% mutate(across(contains('mean'), as.numeric)) %>% mutate(across(contains('LFC'), as.numeric)) %>% mutate(across(contains('p_'), as.numeric))
  return(marker)
}

annotation <- function(data, dataset){
  plan('sequential')

  th0.m.marker <- read.marker('../data/41467_2020_15543_MOESM6_ESM.xlsx', sheet = 3)
  th0.m.marker.up <- th0.m.marker[th0.m.marker$LFC > 1,]
  th0.m.marker.down <- th0.m.marker[th0.m.marker$LFC < -1,]

  th0.n.marker <- read.marker('../data/41467_2020_15543_MOESM6_ESM.xlsx', sheet = 1)
  th0.n.marker.up <- th0.n.marker[th0.n.marker$LFC > 1,]
  th0.n.marker.down <- th0.n.marker[th0.n.marker$LFC < -1,]

  th0.m.marker.5d <- read.marker('../data/41467_2020_15543_MOESM6_ESM.xlsx', sheet = 4)
  th0.m.marker.5d.up <- th0.m.marker.5d[th0.m.marker.5d$LFC > 1,]
  th0.m.marker.5d.down <- th0.m.marker.5d[th0.m.marker.5d$LFC < -1,]

  th0.n.marker.5d <- read.marker('../data/41467_2020_15543_MOESM6_ESM.xlsx', sheet = 2)
  th0.n.marker.5d.up <- th0.n.marker.5d[th0.n.marker.5d$LFC > 1,]
  th0.n.marker.5d.down <- th0.n.marker.5d[th0.n.marker.5d$LFC < -1,]

  th0.marker.up <- intersect(intersect(th0.m.marker.up$gene_name, th0.n.marker.up$gene_name), intersect(th0.m.marker.5d.up$gene_name, th0.n.marker.5d.up$gene_name))
  th0.marker.down <- intersect(intersect(th0.m.marker.down$gene_name, th0.n.marker.down$gene_name), intersect(th0.m.marker.5d.down$gene_name, th0.n.marker.5d.down$gene_name))

  DefaultAssay(data) <- 'SCT'
  data <- AddModuleScore(
    object = data,
    features = list(th0.marker.up),
    name = 'Activated'
  )
  data <- AddModuleScore(
    object = data,
    features = list(th0.marker.down),
    name = 'Resting'
  )
  # p1 <- FeaturePlot(data, 'Activated1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle('Activated')
  # p2 <- FeaturePlot(data, 'Resting1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle('Resting')
  # pdf(paste0('../plots/', dataset, '/activated_resting.pdf'), width = 8, height = 4)
  # print(p1 | p2)
  # dev.off()
  #
  # pdf(paste0('../plots/', dataset, '/naive_memory.pdf'), width = 8, height = 8)
  # print(wrap_plots(check_marker('adt', 'CD4'), check_marker('adt', 'CD8'),
  #                  check_marker('adt', 'CD45RA'), check_marker('adt', 'CD45RO'), ncol = 2))
  # dev.off()
  #
  # pdf(paste0('../plots/', dataset, '/CD4_naive.pdf'), width = 8, height = 8)
  # print(wrap_plots(check_marker('adt', 'CD4'), check_marker('adt', 'CD45RA'),
  #                  check_marker('adt', 'CD62L'), check_marker('rna', 'CCR7'), ncol = 2))
  # dev.off()
  #
  # pdf(paste0('../plots/', dataset, '/CD8_naive.pdf'), width = 8, height = 8)
  # print(wrap_plots(check_marker('adt', 'CD8'), check_marker('adt', 'CD45RA'),
  #                  check_marker('adt', 'CD62L'), check_marker('rna', 'CCR7'), ncol = 2))
  # dev.off()
  #
  # pdf(paste0('../plots/', dataset, '/activated_Th.pdf'), width = 16, height = 8)
  # print(wrap_plots(check_marker('adt', 'CD4'), check_marker('adt', 'CD44'),
  #                  check_marker('adt', 'CD25'), check_marker('adt', 'CD69'),
  #                  check_marker('adt', 'CD71'),
  #                  check_marker('rna', 'SLC3A2'), check_marker('rna', 'SLC7A5'), ncol = 4))
  # dev.off()
  #
  # pdf(paste0('../plots/', dataset, '/Treg.pdf'), width = 16, height = 16)
  # print(wrap_plots(check_marker('adt', 'CD152'),
  #                  check_marker('adt', 'CD279'), check_marker('adt', 'CD39'),
  #                  check_marker('adt', 'CD103'),
  #                  check_marker('adt', 'CD25'), check_marker('adt', 'CD134'),
  #                  check_marker('adt', 'CD137'), check_marker('adt', 'CD223'),
  #                  check_marker('rna', 'CCR7'), check_marker('motif', 'MA0850.1', 'FOXP3'),
  #                  check_marker('rna', 'FOXP3'), check_marker('rna', 'IKZF2'),
  #                  check_marker('rna', 'SMAD3'), check_marker('rna', 'AHR'),
  #                  check_marker('rna', 'CTLA4'), check_marker('rna', 'IL2RA'), ncol = 4))
  # dev.off()
  #
  # pdf(paste0('../plots/', dataset, '/Th17.pdf'), width = 20, height = 16)
  # print(wrap_plots(check_marker('adt', 'CD161'), check_marker('adt', 'CD278'),
  #                  check_marker('adt', 'CD194'), check_marker('adt', 'CD196'),
  #                  check_marker('rna', 'IL21R'), check_marker('rna', 'IL23R'),
  #                  check_marker('rna', 'IL17F'),
  #                  check_marker('rna', 'TNF'), check_marker('rna', 'CCL20'),
  #                  check_marker('motif', 'MA1151.1', 'RORC'),
  #                  check_marker('motif', 'MA0071.1', 'RORA'), check_marker('motif', 'MA0072.1', 'RORA'), check_marker('rna', 'RORA'),
  #                  check_marker('motif', 'MA1634.1', 'BATF'), check_marker('rna', 'BATF'),
  #                  check_marker('motif', 'MA1419.1', 'IRF4'), check_marker('rna', 'IRF4'),
  #                  check_marker('motif', 'MA1419.1', 'IRF4'), check_marker('rna', 'IRF4'), ncol = 5))
  # dev.off()
  #
  # pdf(paste0('../plots/', dataset, '/Th1.pdf'), width = 16, height = 12)
  # print(wrap_plots(check_marker('adt', 'CD183'), check_marker('adt', 'CD195'),
  #                  check_marker('adt', 'CD26'), check_marker('adt', 'CD94'),
  #                  check_marker('adt', 'CD278'),
  #                  check_marker('rna', 'IL18R1'), check_marker('rna', 'IFNG'),
  #                  check_marker('rna', 'LTA'), check_marker('rna', 'TNF'),
  #                  check_marker('motif', 'MA0690.1', 'TBX21'), check_marker('rna', 'TBX21'), ncol = 4))
  # dev.off()
  #
  # pdf(paste0('../plots/', dataset, '/Tfh.pdf'), width = 20, height = 16)
  # print(wrap_plots(check_marker('adt', 'CD185'), check_marker('adt', 'CD84'),
  #                  check_marker('adt', 'CD150'), check_marker('adt', 'CD200'),
  #                  check_marker('adt', 'CD272'), check_marker('adt', 'CD278'),
  #                  check_marker('adt', 'CD279'), check_marker('adt', 'TIGIT'),
  #                  check_marker('adt', 'CD57'), check_marker('adt', 'CD154'),
  #                  check_marker('adt', 'CD304'),
  #                  check_marker('rna', 'IL6R'), check_marker('rna', 'TNFSF8'),
  #                  check_marker('rna', 'IL21R'), check_marker('rna', 'CD40LG'),
  #                  check_marker('motif', 'MA0463.2', 'BCL6'), check_marker('rna', 'BCL6'), ncol = 5))
  # dev.off()
  #
  # pdf(paste0('../plots/', dataset, '/Mono.pdf'), width = 16, height = 12)
  # print(wrap_plots(check_marker('adt', 'CD14'), check_marker('adt', 'CD16'),
  #                  check_marker('adt', 'CD192'), check_marker('adt', 'CD83'),
  #                  check_marker('adt', 'CD123'), check_marker('adt', 'CD99'),
  #                  check_marker('adt', 'CD64'), check_marker('adt', 'CD35'),
  #                  check_marker('adt', 'CD11b'), check_marker('adt', 'CD88'),
  #                  check_marker('adt', 'CX3CR1'),
  #                  check_marker('rna', 'CD83'), ncol = 4))
  # dev.off()
}

#annotation(data, 'downsample')

data$tissue <- 'In-house'

dogma <- readRDS('~/RWorkSpace/GEO_data/DOGMA-seq/Analysis/data_peak_calling_tcell.RDS')
dogma$tissue <- 'Public'

data <- merge(data, dogma)
rm(dogma)

plan('multisession')
options(future.globals.maxSize= 36*1024^3)
DefaultAssay(data) <- 'RNA'
data <- SCTransform(data, method = "glmGamPoi", vars.to.regress = "percent.mt") %>% RunPCA(reduction.name = "pca") %>%
  RunHarmony(group.by.vars = c("tissue"), max.iter.harmony = 20, reduction = 'pca', assay.use = 'SCT',project.dim = FALSE,  reduction.save = "harmony_RNA")
data <- RunUMAP(data, reduction = "harmony_RNA", dims = 1:30, assay = 'SCT', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
data <- FindNeighbors(data, dims = 1:30, reduction = 'harmony_RNA')
data <- FindClusters(data, resolution = 0.2, graph.name = 'SCT_snn', algorithm = 3)
data$merged_rna_clusters <- Idents(data)

DefaultAssay(data) <- 'ADT'
# # we will use all ADT features for dimensional reduction
# # we set a dimensional reduction name to avoid overwriting the
VariableFeatures(data) <- rownames(data[["ADT"]])
data <- NormalizeData(data, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca') %>%
  RunHarmony(group.by.vars = c("tissue"), max.iter.harmony = 20, reduction = 'apca', assay.use = 'ADT',project.dim = FALSE,  reduction.save = "harmony_ADT")
data <- RunUMAP(data, reduction = "harmony_ADT", dims = 1:30, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
data <- FindNeighbors(data, dims = 1:30, reduction = 'harmony_ADT')
data <- FindClusters(data, resolution = 0.2, graph.name = 'ADT_snn', algorithm = 3)
data$merged_adt_clusters <- Idents(data)

DefaultAssay(data) <- "peaks"
data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = 'q75')
data <- RunSVD(data)
data <- RunHarmony(data, group.by.vars = c("tissue"), max.iter.harmony = 20, reduction = 'lsi', assay.use = 'peaks',project.dim = FALSE,  reduction.save = "harmony_peaks")
data <- RunUMAP(data, reduction = "harmony_peaks", dims = 2:30, assay = 'peaks', reduction.name = 'atac.umap', reduction.key = 'atacUMAP_')
data <- FindNeighbors(data, dims = 2:30, reduction = 'harmony_peaks')
data <- FindClusters(data, resolution = 0.2, graph.name = 'peaks_snn', algorithm = 3)
data$merged_atac_clusters <- Idents(data)

# Now run multimodal neighbors and embedding
data <- FindMultiModalNeighbors(object = data,
                                reduction.list = list("harmony_RNA", "harmony_peaks", "harmony_ADT"),
                                dims.list = list(1:30, 2:30, 1:30))
data <- RunUMAP(data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_" )
data <- FindClusters(data, graph.name = "wsnn", resolution = 0.2, algorithm = 3)

data$merged_wnn_clusters <- Idents(data)

data <- FindMultiModalNeighbors(object = data,
                                reduction.list = list("harmony_peaks", "harmony_ADT"),
                                dims.list = list(2:30, 1:30),
                                knn.graph.name = "wknn2",
                                snn.graph.name = "wsnn2",
                                weighted.nn.name = "weighted.nn2")
data <- RunUMAP(data, nn.name = "weighted.nn2", reduction.name = "wnn2.umap", reduction.key = "wnn2UMAP_" )
data <- FindClusters(data, graph.name = "wsnn2", algorithm = 3, resolution = 0.2)
data$merged_wnn2_clusters <- Idents(data)

data <- FindMultiModalNeighbors(object = data,
                                reduction.list = list("harmony_peaks", "harmony_RNA"),
                                dims.list = list(2:30, 1:30),
                                knn.graph.name = "wknn3",
                                snn.graph.name = "wsnn3",
                                weighted.nn.name = "weighted.nn3")
data <- RunUMAP(data, nn.name = "weighted.nn3", reduction.name = "wnn3.umap", reduction.key = "wnn3UMAP_" )
data <- FindClusters(data, graph.name = "wsnn3", algorithm = 3, resolution = 0.2)
data$merged_wnn3_clusters <- Idents(data)

saveRDS(data, file = '../output/downsample/data_downsample_merge_pbmc_tcell.RDS')

pdf('../plots/downsample/merge_pbmc/UMAP_split.pdf', width = 40, height = 32)
DimPlot(data, reduction = 'rna.umap', group.by = 'merged_rna_clusters', split.by = 'sample', ncol = 5, raster = T)
dev.off()

pdf('../plots/downsample/merge_pbmc/UMAP_split_condition.pdf', width = 32, height = 16)
DimPlot(data, reduction = 'rna.umap', group.by = 'merged_rna_clusters', split.by = 'condition', ncol = 4, raster = T)
dev.off()

data$condition <- factor(data$condition, levels = c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2', 'Act_IL1B_IL23_TGFB', 'Act_IL1B_IL23_PGE2_TGFB', 'STIM', 'CTRL'))
p1 <- DimPlot(data[,!data$condition %in% c('CTRL', 'STIM')], reduction = 'rna.umap', split.by = 'condition', ncol = 4, raster = T) + NoLegend()
p2 <- DimPlot(data[,sample(colnames(data)[data$condition %in% c('CTRL', 'STIM')], 10000)], reduction = 'rna.umap', split.by = 'condition', ncol = 2, raster = T) + NoLegend()
p3 <- FeaturePlot(data, 'Activated1', cells = colnames(data)[!data$condition %in% c('CTRL', 'STIM')], reduction = 'rna.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + NoLegend() + ggtitle('Activated')
p4 <- FeaturePlot(data, 'Resting1', cells = colnames(data)[!data$condition %in% c('CTRL', 'STIM')], reduction = 'rna.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + NoLegend() + ggtitle('Resting')
pdf('../plots/downsample/merge_pbmc/UMAP.pdf', width = 32, height = 16)
p1/((p3 | p4) - p2)
dev.off()

data$celltype_updated <- factor(data$celltype_updated, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
                                                  'CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
                                                  'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1',
                                                  'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17',
                                                  'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh',
                                                  'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other',
                                                  'CD8+ Naive (Resting)', 'CD8+ Naive (Activated)',
                                                  'CD8+ Regulatory',
                                                  'CD8+ Memory (Resting)', 'CD8+ Memory (Activated)',
                                                  'MAITs (Resting)', 'MAITs (Activated)', 'Gamma Delta'
))

Idents(data)='celltype_updated'
p5=DimPlot(data, reduction = 'rna.umap', raster = T, cells = colnames(data)[!data$condition %in% c('CTRL', 'STIM')]) + guides(col = guide_legend(ncol = 2, override.aes = list(size=4))) + ggtitle('In-house') + theme(plot.title = element_text(hjust = 0.5))
Idents(data)='tissue'
p6=DimPlot(data, reduction = 'rna.umap', raster = T) + theme(legend.position = 'left')

pdf('../plots/downsample/merge_pbmc/UMAP_updated.pdf', width = 16, height = 12)
(p5 | p6)/p1/((p3 | p4) - p2)
dev.off()

pdf('../plots/downsample/merge_pbmc/UMAP_WNN_split.pdf', width = 40, height = 32)
DimPlot(data, reduction = 'wnn3.umap', group.by = 'merged_wnn3_clusters', split.by = 'sample', ncol = 5, raster = T)
dev.off()

pdf('../plots/downsample/merge_pbmc/UMAP_WNN_split_condition.pdf', width = 32, height = 16)
DimPlot(data, reduction = 'wnn3.umap', group.by = 'merged_wnn3_clusters', split.by = 'condition', ncol = 4, raster = T)
dev.off()

Idents(data)='merged_wnn3_clusters'
data$condition <- factor(data$condition, levels = c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2', 'Act_IL1B_IL23_TGFB', 'Act_IL1B_IL23_PGE2_TGFB', 'STIM', 'CTRL'))
p1 <- DimPlot(data[,!data$condition %in% c('CTRL', 'STIM')], reduction = 'wnn3.umap', split.by = 'condition', ncol = 4, raster = T) + NoLegend()
p2 <- DimPlot(data[,sample(colnames(data)[data$condition %in% c('CTRL', 'STIM')], 10000)], reduction = 'wnn3.umap', split.by = 'condition', ncol = 2, raster = T) + NoLegend()
p3 <- FeaturePlot(data, 'Activated1', cells = colnames(data)[!data$condition %in% c('CTRL', 'STIM')], reduction = 'wnn3.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + NoLegend() + ggtitle('Activated')
p4 <- FeaturePlot(data, 'Resting1', cells = colnames(data)[!data$condition %in% c('CTRL', 'STIM')], reduction = 'wnn3.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + NoLegend() + ggtitle('Resting')
pdf('../plots/downsample/merge_pbmc/UMAP_WNN.pdf', width = 32, height = 16)
p1/((p3 | p4) - p2)
dev.off()

Idents(data)='celltype_updated'
p5=DimPlot(data, reduction = 'wnn3.umap', raster = T, cells = colnames(data)[!data$condition %in% c('CTRL', 'STIM')]) + guides(col = guide_legend(ncol = 2, override.aes = list(size=4))) + ggtitle('In-house') + theme(plot.title = element_text(hjust = 0.5))
Idents(data)='tissue'
p6=DimPlot(data, reduction = 'wnn3.umap', raster = T) + theme(legend.position = 'left')

pdf('../plots/downsample/merge_pbmc/UMAP_WNN_updated.pdf', width = 16, height = 12)
(p5 | p6)/p1/((p3 | p4) - p2)
dev.off()

Idents(data)='merged_wnn3_clusters'
data$condition <- factor(data$condition, levels = c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2', 'Act_IL1B_IL23_TGFB', 'Act_IL1B_IL23_PGE2_TGFB', 'STIM', 'CTRL'))
p1 <- DimPlot(data[,!data$condition %in% c('CTRL', 'STIM')], reduction = 'atac.umap', split.by = 'condition', ncol = 4, raster = T) + NoLegend()
p2 <- DimPlot(data[,sample(colnames(data)[data$condition %in% c('CTRL', 'STIM')], 10000)], reduction = 'atac.umap', split.by = 'condition', ncol = 2, raster = T) + NoLegend()
p3 <- FeaturePlot(data, 'Activated1', cells = colnames(data)[!data$condition %in% c('CTRL', 'STIM')], reduction = 'atac.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + NoLegend() + ggtitle('Activated')
p4 <- FeaturePlot(data, 'Resting1', cells = colnames(data)[!data$condition %in% c('CTRL', 'STIM')], reduction = 'atac.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + NoLegend() + ggtitle('Resting')
pdf('../plots/downsample/merge_pbmc/UMAP_ATAC.pdf', width = 32, height = 16)
p1/((p3 | p4) - p2)
dev.off()

Idents(data)='celltype_updated'
p5=DimPlot(data, reduction = 'atac.umap', raster = T, cells = colnames(data)[!data$condition %in% c('CTRL', 'STIM')]) + guides(col = guide_legend(ncol = 2, override.aes = list(size=4))) + ggtitle('In-house') + theme(plot.title = element_text(hjust = 0.5))
Idents(data)='tissue'
p6=DimPlot(data, reduction = 'atac.umap', raster = T) + theme(legend.position = 'left')

pdf('../plots/downsample/merge_pbmc/UMAP_ATAC_updated.pdf', width = 16, height = 12)
(p5 | p6)/p1/((p3 | p4) - p2)
dev.off()

Idents(data)='tissue'
p6=DimPlot(data, reduction = 'atac.umap', raster = T, cells = colnames(data)[data$tissue %in% c('Public')]) + NoLegend() + ggtitle('Public') + theme(plot.title = element_text(hjust = 0.5))
pdf('../plots/downsample/merge_pbmc/UMAP_ATAC_updated_updated.pdf', width = 16, height = 12)
(p5 | p6)/p1/((p3 | p4) - p2)
dev.off()

Idents(data)='celltype_updated'
p1=DimPlot(data, reduction = 'atac.umap', raster = T, cells = colnames(data)[!data$condition %in% c('CTRL', 'STIM')]) + guides(col = guide_legend(ncol = 2, override.aes = list(size=4))) + ggtitle('In-house') + theme(plot.title = element_text(hjust = 0.5))
Idents(data)='tissue'
p2=DimPlot(data, reduction = 'atac.umap', raster = T, cells = colnames(data)[data$tissue %in% c('Public')]) + NoLegend() + ggtitle('Public') + theme(plot.title = element_text(hjust = 0.5))
Idents(data)='merged_wnn3_clusters'
p3 <- FeaturePlot(data, 'Activated1', cells = colnames(data)[!data$condition %in% c('CTRL', 'STIM')], reduction = 'atac.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + NoLegend() + ggtitle('Activated')
p4 <- FeaturePlot(data, 'Resting1', cells = colnames(data)[!data$condition %in% c('CTRL', 'STIM')], reduction = 'atac.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + NoLegend() + ggtitle('Resting')
p5=p3 + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p3$data[p3$data$Activated1 > median(p3$data$Activated1),], bins = 5)
p6=p4 + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p4$data[p4$data$Resting1 > median(p4$data$Resting1),], bins = 5)
Idents(data)='merged_wnn3_clusters'
p7=DimPlot(data, reduction = 'atac.umap', raster = T, cells = colnames(data)[data$condition %in% c('Act_IL1B_IL23')]) + NoLegend() + ggtitle('Act_IL1B_IL23') + theme(plot.title = element_text(hjust = 0.5)) + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p3$data[p3$data$Activated1 > median(p3$data$Activated1),], bins = 5)
p8=DimPlot(data, reduction = 'atac.umap', raster = T, cells = colnames(data)[data$condition %in% c('Act_IL1B_IL23_PGE2')]) + NoLegend() + ggtitle('Act_IL1B_IL23_PGE2') + theme(plot.title = element_text(hjust = 0.5)) + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p4$data[p4$data$Resting1 > median(p4$data$Resting1),], bins = 5)
p9=DimPlot(data, reduction = 'atac.umap', raster = T, cells = colnames(data)[data$condition %in% c('Act_IL1B_IL23_TGFB')]) + NoLegend() + ggtitle('Act_IL1B_IL23_TGFB') + theme(plot.title = element_text(hjust = 0.5)) + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p3$data[p3$data$Activated1 > median(p3$data$Activated1),], bins = 5)
p10=DimPlot(data, reduction = 'atac.umap', raster = T, cells = colnames(data)[data$condition %in% c('Act_IL1B_IL23_PGE2_TGFB')]) + NoLegend() + ggtitle('Act_IL1B_IL23_PGE2_TGFB') + theme(plot.title = element_text(hjust = 0.5)) + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p4$data[p4$data$Resting1 > median(p4$data$Resting1),], bins = 5)
p11=DimPlot(data, reduction = 'atac.umap', raster = T, cells = colnames(data)[data$condition %in% c('STIM')]) + NoLegend() + ggtitle('STIM') + theme(plot.title = element_text(hjust = 0.5)) + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p3$data[p3$data$Activated1 > median(p3$data$Activated1),], bins = 5)
p12=DimPlot(data, reduction = 'atac.umap', raster = T, cells = colnames(data)[data$condition %in% c('CTRL')]) + NoLegend() + ggtitle('CTRL') + theme(plot.title = element_text(hjust = 0.5)) + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p4$data[p4$data$Resting1 > median(p4$data$Resting1),], bins = 5)


pdf('../plots/downsample/merge_pbmc/UMAP_ATAC_updated_updated_lined.pdf', width = 16, height = 12)
(p1|p2)/(p5|p6|p11|p12)/(p7|p8|p9|p10)
dev.off()

p7=DimPlot(data, reduction = 'atac.umap', raster = T, cells = sample(colnames(data)[data$condition %in% c('Act_IL1B_IL23')], 2000)) + NoLegend() + ggtitle('Act_IL1B_IL23') + theme(plot.title = element_text(hjust = 0.5)) + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p3$data[p3$data$Activated1 > median(p3$data$Activated1),], bins = 5)
p8=DimPlot(data, reduction = 'atac.umap', raster = T, cells = sample(colnames(data)[data$condition %in% c('Act_IL1B_IL23_PGE2')], 2000)) + NoLegend() + ggtitle('Act_IL1B_IL23_PGE2') + theme(plot.title = element_text(hjust = 0.5)) + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p4$data[p4$data$Resting1 > median(p4$data$Resting1),], bins = 5)
p9=DimPlot(data, reduction = 'atac.umap', raster = T, cells = sample(colnames(data)[data$condition %in% c('Act_IL1B_IL23_TGFB')], 2000)) + NoLegend() + ggtitle('Act_IL1B_IL23_TGFB') + theme(plot.title = element_text(hjust = 0.5)) + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p3$data[p3$data$Activated1 > median(p3$data$Activated1),], bins = 5)
p10=DimPlot(data, reduction = 'atac.umap', raster = T, cells = sample(colnames(data)[data$condition %in% c('Act_IL1B_IL23_PGE2_TGFB')], 2000)) + NoLegend() + ggtitle('Act_IL1B_IL23_PGE2_TGFB') + theme(plot.title = element_text(hjust = 0.5)) + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p4$data[p4$data$Resting1 > median(p4$data$Resting1),], bins = 5)
p11=DimPlot(data, reduction = 'atac.umap', raster = T, cells = sample(colnames(data)[data$condition %in% c('STIM')], 2000)) + NoLegend() + ggtitle('STIM') + theme(plot.title = element_text(hjust = 0.5)) + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p3$data[p3$data$Activated1 > median(p3$data$Activated1),], bins = 5)
p12=DimPlot(data, reduction = 'atac.umap', raster = T, cells = sample(colnames(data)[data$condition %in% c('CTRL')], 2000)) + NoLegend() + ggtitle('CTRL') + theme(plot.title = element_text(hjust = 0.5)) + geom_density_2d(aes(x = atacUMAP_1, y = atacUMAP_2), data = p4$data[p4$data$Resting1 > median(p4$data$Resting1),], bins = 5)


pdf('../plots/downsample/merge_pbmc/UMAP_ATAC_updated_updated_lined_downsampled.pdf', width = 16, height = 12)
(p1|p2)/(p5|p6|p11|p12)/(p7|p8|p9|p10)
dev.off()
