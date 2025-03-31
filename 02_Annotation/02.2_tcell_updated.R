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

# load('../output/data_harmony_qc_tcell.RData')
# 
# data <- data[,data$harmony_atac_clusters != 8]
# 
# plot <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = 'predicted.celltype.l1')
# cell1 <- CellSelector(plot)
# cell2 <- CellSelector(plot)
# DimPlot(data, reduction = 'rna.umap', cells.highlight = c(cell1, cell2))
# DimPlot(data, reduction = 'adt.umap', cells.highlight = c(cell1, cell2))
# DimPlot(data, reduction = 'atac.umap', cells.highlight = c(cell1, cell2))
# DimPlot(data, reduction = 'wnn.umap', cells.highlight = c(cell1, cell2))
# 
# data <- data[,!colnames(data) %in% c(cell1, cell2)]
# DimPlot(data, reduction = 'atac.umap', cells.highlight = c(cell1, cell2))
# 
# save(data, file = '../output/data_harmony_qc_tcell_updated.RData')

load('../output/data_harmony_qc_tcell_updated.RData')

plan('multicore')
options(future.globals.maxSize= 60*1024^3)

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

save(data, file = '../output/data_harmony_qc_tcell_updated.RData')

dataset <- 'qc_tcell_updated'

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