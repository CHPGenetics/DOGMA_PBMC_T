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

data <- readRDS('../output/tcell_annotated.RDS')

p1 <- FeaturePlot(data, 'Activated1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle('Activated')
p2 <- FeaturePlot(data, 'Resting1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle('Resting') 
pdf('../plots/annotation_all/activated_resting.pdf', width = 16, height = 8)
print(p1 | p2)
dev.off()

data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = c(seq(0.02, 0.2, 0.02)))

pdf(paste0('../plots/annotation_all/tree_WNN.pdf'), width = 10, height = 10)
print(clustree(data@meta.data, prefix = 'wsnn_res.'))
dev.off()

data <- FindClusters(data, graph.name = "SCT_snn", algorithm = 3, resolution = c(seq(0.02, 0.2, 0.02)))

pdf(paste0('../plots/annotation_all/tree_RNA.pdf'), width = 10, height = 10)
print(clustree(data@meta.data, prefix = 'SCT_snn_res.'))
dev.off()

data <- FindClusters(data, graph.name = "ADT_snn", algorithm = 3, resolution = c(seq(0.02, 0.2, 0.02)))

pdf(paste0('../plots/annotation_all/tree_ADT.pdf'), width = 10, height = 10)
print(clustree(data@meta.data, prefix = 'ADT_snn_res.'))
dev.off()

data <- FindClusters(data, graph.name = "peaks_snn", algorithm = 3, resolution = c(seq(0.02, 0.2, 0.02)))

pdf(paste0('../plots/annotation_all/tree_ATAC.pdf'), width = 10, height = 10)
print(clustree(data@meta.data, prefix = 'peaks_snn_res.'))
dev.off()

data <- FindClusters(data, graph.name = "wsnn2", algorithm = 3, resolution = c(seq(0.02, 0.2, 0.02)))

pdf(paste0('../plots/annotation_all/tree_WNN2.pdf'), width = 10, height = 10)
print(clustree(data@meta.data, prefix = 'wsnn2_res.'))
dev.off()

saveRDS(data, file = '../output/tcell_annotated_updated.RDS')