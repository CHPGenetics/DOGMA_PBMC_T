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

load('../output/tcell_chrom.RData')

# data_list <- SplitObject(data, split.by = "orig.ident")
# 
# library(parallel)
# data_list_analyzed <- mclapply(data_list, FUN = function(data){
#   dataset <- names(table(data$orig.ident))[1]
#   print(dataset)
#   dir.create(file.path(paste0('../plots/annotation_all'), dataset), showWarnings = FALSE)
#   
#   DefaultAssay(data) <- "RNA"
#   data <- SCTransform(data, method = "glmGamPoi", vars.to.regress = "percent.mt") %>% RunPCA(reduction.name = "pca")
#   data <- RunUMAP(data, reduction = "pca", dims = 1:30, assay = 'SCT', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
#   data <- FindNeighbors(data, dims = 1:30, reduction = 'pca')
#   data <- FindClusters(data, resolution = 0.2, graph.name = 'SCT_snn', algorithm = 3)
#   data$rna_clusters <- Idents(data)
#   
#   DefaultAssay(data) <- 'ADT'
#   # we will use all ADT features for dimensional reduction
#   # we set a dimensional reduction name to avoid overwriting the
#   VariableFeatures(data) <- rownames(data[["ADT"]])
#   data <- NormalizeData(data, normalization.method = 'CLR', margin = 2) %>%
#     ScaleData() %>% RunPCA(reduction.name = 'apca')
#   data <- RunUMAP(data, reduction = "apca", dims = 1:30, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
#   data <- FindNeighbors(data, dims = 1:30, reduction = 'apca')
#   data <- FindClusters(data, resolution = 0.2, graph.name = 'ADT_snn', algorithm = 3)
#   data$adt_clusters <- Idents(data)
#   
#   DefaultAssay(data) <- "peaks"
#   data <- RunTFIDF(data)
#   data <- FindTopFeatures(data, min.cutoff = 'q0')
#   data <- RunSVD(data)
#   data <- RunUMAP(data, reduction = "lsi", dims = 2:30, assay = 'peaks', reduction.name = 'atac.umap', reduction.key = 'atacUMAP_')
#   data <- FindNeighbors(data, dims = 2:30, reduction = 'lsi')
#   data <- FindClusters(data, resolution = 0.2, graph.name = 'peaks_snn', algorithm = 3)
#   data$atac_clusters <- Idents(data)
#   
#   # Now run multimodal neighbors and embedding
#   data <- FindMultiModalNeighbors(object = data,
#                                   reduction.list = list("pca", "lsi", "apca"),
#                                   dims.list = list(1:30, 2:30, 1:30))
#   data <- RunUMAP(data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_" )
#   data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = 0.2)
#   
#   data$wnn_clusters <- Idents(data)
#   
#   data <- FindMultiModalNeighbors(object = data,
#                                   reduction.list = list("lsi", "apca"),
#                                   dims.list = list(2:30, 1:30),
#                                   knn.graph.name = "wknn2",
#                                   snn.graph.name = "wsnn2",
#                                   weighted.nn.name = "weighted.nn2")
#   data <- RunUMAP(data, nn.name = "weighted.nn2", reduction.name = "wnn2.umap", reduction.key = "wnn2UMAP_" )
#   data <- FindClusters(data, graph.name = "wsnn2", algorithm = 3, resolution = 0.2)
#   data$wnn2_clusters <- Idents(data)
#   
#   reference <- readRDS('../output/tcell_downsample_annotated.RDS')
#   
#   DefaultAssay(data) <- 'SCT'
#   DefaultAssay(reference) <- 'SCT'
#   
#   anchors.rna <- FindTransferAnchors(
#     reference = reference,
#     query = data,
#     normalization.method = "SCT",
#     reference.reduction = "rna.spca",
#     reduction = 'pcaproject',
#     dims = 1:30,
#     reference.assay = 'SCT',
#     query.assay = 'SCT'
#   )
#   
#   anchors.adt <- FindTransferAnchors(
#     reference = reference,
#     query = data,
#     reference.reduction = "adt.spca",
#     reduction = 'pcaproject',
#     dims = 1:30,
#     reference.assay = 'ADT',
#     query.assay = 'ADT'
#   )
#   
#   anchors.atac <- FindTransferAnchors(
#     reference = reference,
#     query = data,
#     reference.reduction = "lsi",
#     reduction = 'lsiproject',
#     dims = 2:30,
#     reference.assay = 'peaks',
#     query.assay = 'peaks'
#   )
#   
#   gc()
#   
#   predictions.rna <- TransferData(anchorset = anchors.rna, refdata = reference$celltype, weight.reduction = "pcaproject")
#   predictions.adt <- TransferData(anchorset = anchors.adt, refdata = reference$celltype, weight.reduction = "pcaproject")
#   predictions.atac <- TransferData(anchorset = anchors.atac, refdata = reference$celltype, weight.reduction = "lsiproject")
#   
#   gc()
#   
#   predictions.max <- data.table(rna = rowMax(as.matrix(predictions.rna[,-1])),
#                                 adt = rowMax(as.matrix(predictions.adt[,-1])),
#                                 atac = rowMax(as.matrix(predictions.atac[,-1])))
#   rownames(predictions.max) <- rownames(predictions.rna)
#   predictions.max.label <- apply(predictions.max, 1, which.max)
#   names(predictions.max.label) <- rownames(predictions.rna)
#   
#   predictions <- rbind(predictions.rna[names(predictions.max.label[predictions.max.label == 1]),],
#                        predictions.adt[names(predictions.max.label[predictions.max.label == 2]),],
#                        predictions.atac[names(predictions.max.label[predictions.max.label == 3]),])
#   predictions <- predictions[rownames(predictions.rna),1]
#   
#   data$celltype <- factor(predictions, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
#                                                                       'CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
#                                                                       'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1', 
#                                                                       'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17', 
#                                                                       'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh', 
#                                                                       'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other',
#                                                                       'CD8+ Naive (Resting)', 'CD8+ Naive (Activated)',
#                                                                       'CD8+ Regulatory',
#                                                                       'CD8+ Memory (Resting)', 'CD8+ Memory (Activated)',
#                                                                       'MAITs (Resting)', 'MAITs (Activated)', 'Gamma Delta'
#   ))
#   
#   data$celltype.rna <- factor(predictions.rna$predicted.id, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
#                                                                                            'CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
#                                                                                            'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1', 
#                                                                                            'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17', 
#                                                                                            'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh', 
#                                                                                            'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other',
#                                                                                            'CD8+ Naive (Resting)', 'CD8+ Naive (Activated)',
#                                                                                            'CD8+ Regulatory',
#                                                                                            'CD8+ Memory (Resting)', 'CD8+ Memory (Activated)',
#                                                                                            'MAITs (Resting)', 'MAITs (Activated)', 'Gamma Delta'
#   ))
#   
#   data$celltype.adt <- factor(predictions.adt$predicted.id, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
#                                                                                            'CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
#                                                                                            'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1', 
#                                                                                            'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17', 
#                                                                                            'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh', 
#                                                                                            'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other',
#                                                                                            'CD8+ Naive (Resting)', 'CD8+ Naive (Activated)',
#                                                                                            'CD8+ Regulatory',
#                                                                                            'CD8+ Memory (Resting)', 'CD8+ Memory (Activated)',
#                                                                                            'MAITs (Resting)', 'MAITs (Activated)', 'Gamma Delta'
#   ))
#   
#   data$celltype.atac <- factor(predictions.atac$predicted.id, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
#                                                                                              'CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
#                                                                                              'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1', 
#                                                                                              'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17', 
#                                                                                              'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh', 
#                                                                                              'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other',
#                                                                                              'CD8+ Naive (Resting)', 'CD8+ Naive (Activated)',
#                                                                                              'CD8+ Regulatory',
#                                                                                              'CD8+ Memory (Resting)', 'CD8+ Memory (Activated)',
#                                                                                              'MAITs (Resting)', 'MAITs (Activated)', 'Gamma Delta'
#   ))
#   
#   p1 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = 'wnn_clusters')
#   p2 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "celltype")
#   pdf(paste0('../plots/annotation_all/', dataset, '/annotated_predicted_WNN.pdf'), width = 20, height = 8)
#   print(p1 | p2)
#   dev.off()
#   
#   p1 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = 'rna_clusters')
#   p2 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = "celltype")
#   pdf(paste0('../plots/annotation_all/', dataset, '/annotated_predicted_RNA.pdf'), width = 20, height = 8)
#   print(p1 | p2)
#   dev.off()
#   
#   p1 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = 'adt_clusters')
#   p2 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = "celltype")
#   pdf(paste0('../plots/annotation_all/', dataset, '/annotated_predicted_ADT.pdf'), width = 20, height = 8)
#   print(p1 | p2)
#   dev.off()
#   
#   p1 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = 'atac_clusters')
#   p2 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = "celltype")
#   pdf(paste0('../plots/annotation_all/', dataset, '/annotated_predicted_ATAC.pdf'), width = 20, height = 8)
#   print(p1 | p2)
#   dev.off()
#   
#   p1 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = 'wnn2_clusters')
#   p2 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = "celltype")
#   pdf(paste0('../plots/annotation_all/', dataset, '/annotated_predicted_WNN2.pdf'), width = 20, height = 8)
#   print(p1 | p2)
#   dev.off()
#   
#   write.csv(data@meta.data, file = paste0('../output/annotation_all/', dataset, '.csv'), quote = F, col.names = T, row.names = T)
#   saveRDS(data, file = paste0('../output/annotation_all/', dataset, '.RDS'))
# })

files <- list.files('../output/annotation_all/', '.csv')
meta_list <- lapply(files, FUN = function(file){
  meta <- read.csv(paste0('../output/annotation_all/', file), row.names = 'X')
  return(meta)
})
meta <- Reduce(rbind, meta_list)

meta <- meta[colnames(data),]
data$celltype <- meta$celltype
data$celltype.rna <- meta$celltype.rna
data$celltype.adt <- meta$celltype.adt
data$celltype.atac <- meta$celltype.atac

table(data$celltype)

Idents(data) <- 'celltype'

data$celltype <- factor(data$celltype, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
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

data$celltype.rna <- factor(data$celltype.rna, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
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

data$celltype.adt <- factor(data$celltype.adt, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
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

data$celltype.atac <- factor(data$celltype.atac, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
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

table(data$celltype)
p1 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = 'harmony_wnn_clusters')
p2 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "celltype")
pdf(paste0('../plots/annotation_all/annotated_predicted_WNN.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = 'harmony_rna_clusters')
p2 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = "celltype")
pdf(paste0('../plots/annotation_all/annotated_predicted_RNA.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = 'harmony_adt_clusters')
p2 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = "celltype")
pdf(paste0('../plots/annotation_all/annotated_predicted_ADT.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = 'harmony_atac_clusters')
p2 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = "celltype")
pdf(paste0('../plots/annotation_all/annotated_predicted_ATAC.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = 'harmony_wnn2_clusters')
p2 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = "celltype")
pdf(paste0('../plots/annotation_all/annotated_predicted_WNN2.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

write.csv(data@meta.data, file = '../output/tcell_annotated.csv', quote = F, col.names = T, row.names = T)
saveRDS(data, file = '../output/tcell_annotated.RDS')