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

load('../output/data_macs2.RData')

data$cryopreservation <- 'Fresh'
data$cryopreservation[!data$orig.ident %in% c('Duerr_20210419_DOGMAseq_DIG', 
                                              'Duerr_20210610_DOGMAseq-1', 'Duerr_20210610_DOGMAseq-2',
                                              'Duerr_20210826_DOGMAseq-1', 'Duerr_20210826_DOGMAseq-2',
                                              'Duerr_20210831_DOGMAseq-1', 'Duerr_20210831_DOGMAseq-2')] <- 'Cryopreserved'

Idents(data) <- 'cryopreservation'
DefaultAssay(data) <- 'ADT'

pdf('../plots/macs2_updated/ridge.pdf', width = 20, height = 20)
for (i in 1:6){
  print(RidgePlot(data, features = rownames(data)[(i-1)*25 + 1:25], ncol = 5))
}
print(RidgePlot(data, features = rownames(data)[151:163], ncol = 5))
dev.off()

data@assays$ADT@scale.data <- as.matrix(log10(1 + data@assays$ADT@counts))

pdf('../plots/macs2_updated/ridge_log10.pdf', width = 20, height = 20)
for (i in 1:6){
  print(RidgePlot(data, features = rownames(data)[(i-1)*25 + 1:25], ncol = 5, slot = 'scale.data'))
}
print(RidgePlot(data, features = rownames(data)[151:163], ncol = 5, slot = 'scale.data'))
dev.off()

plan('multicore')
options(future.globals.maxSize= 60*1024^3)

# plan('sequential')
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

save(data, file = '../output/data_macs2_harmony_updated.RData')

pdf('../plots/macs2_updated/UMAP_macs2_harmony_weighted_cluster_split.pdf', width = 40, height = 40)
DimPlot(data, reduction = 'wnn.umap', split.by = 'sample', label = T, repel = TRUE, label.size = 4, ncol = 5)
dev.off()

pdf('../plots/macs2_updated/UMAP_macs2_harmony_weighted_cluster_split_condition.pdf', width = 32, height = 8)
DimPlot(data, reduction = 'wnn.umap', split.by = 'condition', label = T, repel = TRUE, label.size = 4, ncol = 4)
dev.off()

pdf('../plots/macs2_updated/UMAP_macs2_harmony_RNA_cluster_split.pdf', width = 40, height = 40)
DimPlot(data, reduction = 'rna.umap', split.by = 'sample', label = T, repel = TRUE, label.size = 4, ncol = 5)
dev.off()

pdf('../plots/macs2_updated/UMAP_macs2_harmony_RNA_cluster_split_condition.pdf', width = 32, height = 8)
DimPlot(data, reduction = 'rna.umap', split.by = 'condition', label = T, repel = TRUE, label.size = 4, ncol = 4)
dev.off()

pdf('../plots/macs2_updated/UMAP_macs2_harmony_ADT_cluster_split.pdf', width = 40, height = 40)
DimPlot(data, reduction = 'adt.umap', split.by = 'sample', label = T, repel = TRUE, label.size = 4, ncol = 5)
dev.off()

pdf('../plots/macs2_updated/UMAP_macs2_harmony_ADT_cluster_split_condition.pdf', width = 32, height = 8)
DimPlot(data, reduction = 'adt.umap', split.by = 'condition', label = T, repel = TRUE, label.size = 4, ncol = 4)
dev.off()

pdf('../plots/macs2_updated/UMAP_macs2_harmony_ATAC_cluster_split.pdf', width = 40, height = 40)
DimPlot(data, reduction = 'atac.umap', split.by = 'sample', label = T, repel = TRUE, label.size = 4, ncol = 5)
dev.off()

pdf('../plots/macs2_updated/UMAP_macs2_harmony_ATAC_cluster_split_condition.pdf', width = 32, height = 8)
DimPlot(data, reduction = 'atac.umap', split.by = 'condition', label = T, repel = TRUE, label.size = 4, ncol = 4)
dev.off()

pdf('../plots/macs2_updated/UMAP_macs2_harmony_sample_condition.pdf', width = 16, height = 8)
p1 <- DimPlot(data, reduction = 'wnn.umap', group.by = 'sample', label = F, repel = TRUE, label.size = 4)
p2 <- DimPlot(data, reduction = 'wnn.umap', group.by = 'condition', label = F, repel = TRUE, label.size = 4)
p1 | p2
dev.off()

pdf('../plots/macs2_updated/UMAP_macs2_harmony_sample_condition_split.pdf', width = 32, height = 8)
p1 <- DimPlot(data, reduction = 'wnn.umap', group.by = 'sample', split.by = 'condition', label = F, repel = TRUE, label.size = 4, ncol = 4)
p1
dev.off()

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

p1 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE)
p2 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1") 
p3 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l2") 
pdf('../plots/macs2_updated/predicted_macs2.pdf', width = 24, height = 8)
p1 | p2 | p3
dev.off()

p1 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE)
p2 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1") 
pdf('../plots/macs2_updated/predicted_2_macs2.pdf', width = 16, height = 8)
p1 | p2
dev.off()

library(ggplot2)
data_plot <- table(data$predicted.celltype.l1, data$condition)
data_plot <- reshape2::melt(prop.table(data_plot, 2))
colnames(data_plot) <- c('Clusters', 'Conditions', 'Proportions')
data_plot$Clusters <- as.factor(data_plot$Clusters)
data_plot$Proportions <- data_plot$Proportions*100

p1 <- ggplot(data_plot, aes(x = Clusters, y = Proportions, col = Conditions, fill = Conditions))+
  geom_col(position = 'dodge')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0))
pdf('../plots/macs2_updated/proportion_predicted_slide_macs2.pdf', width = 8, height = 4)
p1
dev.off()

data_plot <- table(data$predicted.celltype.l2, data$condition)
data_plot <- reshape2::melt(prop.table(data_plot, 2))
colnames(data_plot) <- c('Clusters', 'Conditions', 'Proportions')
data_plot$Clusters <- as.factor(data_plot$Clusters)
data_plot$Proportions <- data_plot$Proportions*100

p1 <- ggplot(data_plot, aes(x = Clusters, y = Proportions, col = Conditions, fill = Conditions))+
  geom_col(position = 'dodge')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0))
pdf('../plots/macs2_updated/proportion_predicted_l2_slide_macs2.pdf', width = 8, height = 4)
p1
dev.off()

save(data, file = '../output/data_macs2_harmony_updated.RData')