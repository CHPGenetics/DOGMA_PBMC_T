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

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')

load('../output/tcell_downsample.RData')

plan('multicore')
options(future.globals.maxSize= 60*1024^3)

dataset <- 'tcell_downsample'
data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = c(seq(0.02, 0.2, 0.02)))

pdf(paste0('../plots/', dataset, '/tree_WNN.pdf'), width = 10, height = 10)
print(clustree(data@meta.data, prefix = 'wsnn_res.'))
dev.off()

data <- FindClusters(data, graph.name = "SCT_snn", resolution = c(seq(0.02, 0.2, 0.02)))

pdf(paste0('../plots/', dataset, '/tree_RNA.pdf'), width = 10, height = 10)
print(clustree(data@meta.data, prefix = 'SCT_snn_res.'))
dev.off()

data <- FindClusters(data, graph.name = "ADT_snn", resolution = c(seq(0.02, 0.2, 0.02)))

pdf(paste0('../plots/', dataset, '/tree_ADT.pdf'), width = 10, height = 10)
print(clustree(data@meta.data, prefix = 'ADT_snn_res.'))
dev.off()

data <- FindClusters(data, graph.name = "peaks_snn", algorithm = 3, resolution = c(seq(0.02, 0.2, 0.02)))

pdf(paste0('../plots/', dataset, '/tree_ATAC.pdf'), width = 10, height = 10)
print(clustree(data@meta.data, prefix = 'peaks_snn_res.'))
dev.off()

data <- FindClusters(data, graph.name = "wsnn2", algorithm = 3, resolution = c(seq(0.02, 0.2, 0.02)))

pdf(paste0('../plots/', dataset, '/tree_WNN2.pdf'), width = 10, height = 10)
print(clustree(data@meta.data, prefix = 'wsnn2_res.'))
dev.off()

p1 <- FeaturePlot(data, 'percent.mt', reduction = 'wnn.umap')
p2 <- FeaturePlot(data, 'percent.mt', reduction = 'rna.umap')
p3 <- FeaturePlot(data, 'percent.mt', reduction = 'adt.umap')
p4 <- FeaturePlot(data, 'percent.mt', reduction = 'atac.umap')
p5 <- FeaturePlot(data, 'percent.mt', reduction = 'wnn2.umap')
pdf(paste0('../plots/', dataset, '/mito.pdf'), width = 24, height = 16)
print(wrap_plots(p1, p2, p3, p4, p5, ncol = 3))
dev.off()

################
reference <- readRDS('../../Time_series/output/reference_qc.RDS')

DefaultAssay(data) <- 'SCT'
DefaultAssay(reference) <- 'SCT'

anchors.rna <- FindTransferAnchors(
  reference = reference,
  query = data,
  normalization.method = "SCT",
  reference.reduction = "rna.spca",
  reduction = 'pcaproject',
  dims = 1:30,
  reference.assay = 'SCT',
  query.assay = 'SCT'
)

anchors.adt <- FindTransferAnchors(
  reference = reference,
  query = data,
  reference.reduction = "adt.spca",
  reduction = 'pcaproject',
  dims = 1:30,
  reference.assay = 'ADT',
  query.assay = 'ADT'
)

anchors.atac <- FindTransferAnchors(
  reference = reference,
  query = data,
  reference.reduction = "lsi",
  reduction = 'lsiproject',
  dims = 2:30,
  reference.assay = 'peaks',
  query.assay = 'peaks'
)

gc()

predictions.rna <- TransferData(anchorset = anchors.rna, refdata = reference$celltype, weight.reduction = "pcaproject")
predictions.adt <- TransferData(anchorset = anchors.adt, refdata = reference$celltype, weight.reduction = "pcaproject")
predictions.atac <- TransferData(anchorset = anchors.atac, refdata = reference$celltype, weight.reduction = "lsiproject")

gc()

predictions.max <- data.table(rna = rowMax(as.matrix(predictions.rna[,-1])),
                              adt = rowMax(as.matrix(predictions.adt[,-1])),
                              atac = rowMax(as.matrix(predictions.atac[,-1])))
rownames(predictions.max) <- rownames(predictions.rna)
predictions.max.label <- apply(predictions.max, 1, which.max)
names(predictions.max.label) <- rownames(predictions.rna)

predictions <- rbind(predictions.rna[names(predictions.max.label[predictions.max.label == 1]),],
                     predictions.adt[names(predictions.max.label[predictions.max.label == 2]),],
                     predictions.atac[names(predictions.max.label[predictions.max.label == 3]),])
predictions <- predictions[rownames(predictions.rna),1]

data$annotated_predicted.celltype <- factor(predictions, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
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

data$annotated_predicted.celltype.rna <- factor(predictions.rna$predicted.id, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
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

data$annotated_predicted.celltype.adt <- factor(predictions.adt$predicted.id, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
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

data$annotated_predicted.celltype.atac <- factor(predictions.atac$predicted.id, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
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

p1 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = 'wnn_clusters')
p2 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
pdf(paste0('../plots/', dataset, '/annotated_predicted_WNN.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = 'rna_clusters')
p2 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
pdf(paste0('../plots/', dataset, '/annotated_predicted_RNA.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = 'adt_clusters')
p2 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
pdf(paste0('../plots/', dataset, '/annotated_predicted_ADT.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = 'atac_clusters')
p2 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
pdf(paste0('../plots/', dataset, '/annotated_predicted_ATAC.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = 'wnn2_clusters')
p2 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
pdf(paste0('../plots/', dataset, '/annotated_predicted_WNN2.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

save(data, file = '../output/tcell_downsample_predicted.RData')
