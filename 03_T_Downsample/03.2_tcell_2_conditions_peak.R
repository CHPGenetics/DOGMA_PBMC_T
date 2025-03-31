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

data = readRDS('../output/tcell_annotated_updated_2_conditions.RDS')

peak = readRDS('../output/pseudo_bulk/peak.RDS')
peak = peak[,colnames(data)]

DefaultAssay(data) <- 'peaks'
data[['peaks']] <- CreateChromatinAssay(counts = peak,
                                         fragments = Fragments(data),
                                         annotation = Annotation(data))

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

dataset <- '2_conditions_peak'

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

saveRDS(data, file = '../output/tcell_annotated_updated_2_conditions_peak.RDS')

for (i in 1:20){
  celltype <- levels(data$celltype_updated)[i]
  assign(paste0('p',i),
         DimPlot(data, reduction = 'atac.umap', cells.highlight = colnames(data)[data$celltype_updated == celltype]) + NoLegend() + ggtitle(celltype) + theme(plot.title = element_text(hjust = 0.5)))
}

pdf('../plots/2_conditions_peak/celltype.pdf', width = 20, height = 16)
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20, ncol = 5)
dev.off()
