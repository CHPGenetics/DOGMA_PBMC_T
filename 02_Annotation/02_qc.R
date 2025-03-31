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

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')

load('../output/data_macs2_harmony_updated.RData')

plan('multicore')
options(future.globals.maxSize= 60*1024^3)

SumFeatureSet <- function (object, pattern = NULL, features = NULL, col.name = NULL, assay = NULL){
  assay <- assay %||% DefaultAssay(object = object)
  if (!is.null(x = features) && !is.null(x = pattern)) {
    warning("Both pattern and features provided. Pattern is being ignored.")
  }
  features <- features %||% grep(pattern = pattern, x = rownames(x = object[[assay]]),
                                 value = TRUE)
  percent.featureset <- colSums(x = GetAssayData(object = object,
                                                 assay = assay, slot = "counts")[features, , drop = FALSE])
  if (!is.null(x = col.name)) {
    object <- AddMetaData(object = object, metadata = percent.featureset,
                          col.name = col.name)
    return(object)
  }
  return(percent.featureset)
}

data[['sum.ctrl']] <- SumFeatureSet(data, 'Ctrl', assay = 'ADT')

Idents(data) <- 'sample'
p1 <- VlnPlot(data, 'nFeature_RNA', pt.size = 0, log = T) + NoLegend() + ggtitle('') + xlab('') + ylab('# of genes')
p2 <- VlnPlot(data, 'nCount_RNA', pt.size = 0, log = T) + NoLegend() + ggtitle('') + xlab('') + ylab('# of UMIs')
p3 <- VlnPlot(data, 'percent.mt', pt.size = 0) + NoLegend() + ggtitle('') + xlab('') + ylab('%UMIs from mtRNA')
p4 <- VlnPlot(data, 'sum.ctrl', pt.size = 0) + ylim(c(0,100)) + NoLegend() + ggtitle('') + xlab('') + ylab('# of tags from isotype controls')
p5 <- VlnPlot(data, 'atac_fragments', pt.size = 0, log = T) + NoLegend() + ggtitle('') + xlab('') + ylab('# of fragments')
p6 <- VlnPlot(data, 'TSS.enrichment', pt.size = 0, log = T) + NoLegend() + ggtitle('') + xlab('') + ylab('TSS enrichment score')
pdf('../plots/qc/violin_QC_before.pdf', width = 24, height = 4)
print(wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 6))
dev.off()

table(data$percent.mt < 20 & data$nFeature_RNA > 200 & data$nCount_RNA > 1000 & data$TSS.enrichment > 4 & data$atac_fragments > 1000 & data$sum.ctrl < 50)
# FALSE   TRUE 
# 19922 169893

p1 <- DimPlot(data, reduction = 'wnn.umap', cells.highlight = colnames(data)[!(data$percent.mt < 20 & data$nFeature_RNA > 200 & data$nCount_RNA > 1000 & data$TSS.enrichment > 4 & data$atac_fragments > 1000 & data$sum.ctrl < 50)])
p2 <- DimPlot(data, reduction = 'rna.umap', cells.highlight = colnames(data)[!(data$percent.mt < 20 & data$nFeature_RNA > 200 & data$nCount_RNA > 1000 & data$TSS.enrichment > 4 & data$atac_fragments > 1000 & data$sum.ctrl < 50)])
p3 <- DimPlot(data, reduction = 'adt.umap', cells.highlight = colnames(data)[!(data$percent.mt < 20 & data$nFeature_RNA > 200 & data$nCount_RNA > 1000 & data$TSS.enrichment > 4 & data$atac_fragments > 1000 & data$sum.ctrl < 50)])
p4 <- DimPlot(data, reduction = 'atac.umap', cells.highlight = colnames(data)[!(data$percent.mt < 20 & data$nFeature_RNA > 200 & data$nCount_RNA > 1000 & data$TSS.enrichment > 4 & data$atac_fragments > 1000 & data$sum.ctrl < 50)])

pdf('../plots/qc/UMAP_QC.pdf', width = 32, height = 8)
wrap_plots(p1, p2, p3, p4, ncol = 4)
dev.off()

data <- data[,data$percent.mt < 20 & data$nFeature_RNA > 200 & data$nCount_RNA > 1000 & data$TSS.enrichment > 4 & data$atac_fragments > 1000 & data$sum.ctrl < 50]

p1 <- VlnPlot(data, 'nFeature_RNA', pt.size = 0, log = T) + NoLegend() + ggtitle('') + xlab('') + ylab('# of genes')
p2 <- VlnPlot(data, 'nCount_RNA', pt.size = 0, log = T) + NoLegend() + ggtitle('') + xlab('') + ylab('# of UMIs')
p3 <- VlnPlot(data, 'percent.mt', pt.size = 0) + NoLegend() + ggtitle('') + xlab('') + ylab('%UMIs from mtRNA')
p4 <- VlnPlot(data, 'sum.ctrl', pt.size = 0) + NoLegend() + ggtitle('') + xlab('') + ylab('# of tags from isotype controls')
p5 <- VlnPlot(data, 'atac_fragments', pt.size = 0, log = T) + NoLegend() + ggtitle('') + xlab('') + ylab('# of fragments')
p6 <- VlnPlot(data, 'TSS.enrichment', pt.size = 0, log = T) + NoLegend() + ggtitle('') + xlab('') + ylab('TSS enrichment score')
pdf('../plots/qc/violin_QC_after.pdf', width = 24, height = 4)
print(wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 6))
dev.off()

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

reference <- readRDS('~/RWorkSpace/CITE-seq/Seurat/PBMC.RDS')

DefaultAssay(data) <- 'SCT'

anchors <- FindTransferAnchors(
  reference = reference,
  query = data,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

data <- MapQuery(
  anchorset = anchors,
  query = data,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    celltype.l3 = "celltype.l3"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

save(data, file = '../output/data_harmony_qc.RData')

dataset <- 'qc'

p1 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = 'harmony_wnn_clusters')
p2 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1")
pdf(paste0('../plots/', dataset, '/predicted_WNN.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = 'harmony_rna_clusters')
p2 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1")
pdf(paste0('../plots/', dataset, '/predicted_RNA.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = 'harmony_adt_clusters')
p2 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1")
pdf(paste0('../plots/', dataset, '/predicted_ADT.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = 'harmony_atac_clusters')
p2 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1")
pdf(paste0('../plots/', dataset, '/predicted_ATAC.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()