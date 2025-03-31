library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
set.seed(1)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/')

addArchRThreads(threads = 64)

proj <- loadArchRProject('DOGMA_filtered_multiome_2/')

# data <- read.csv('../tcell_annotated_updated.csv', row.names = 'X')
# data <- data[data$condition %in% c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'),]
# 
# table(proj$cellNames %in% gsub('_', '#', rownames(data)))
# proj <- proj[proj$cellNames %in% gsub('_', '#', rownames(data)),]
# 
# proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered_multiome_2')

proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC",
  force = T
)

proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA",
  force = T
)

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "LSI_ATAC",
  name = "Harmony_ATAC",
  groupBy = c("sample"),
  force = T
)

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "LSI_RNA",
  name = "Harmony_RNA",
  groupBy = c("sample"),
  force = T
)
#Combined Dims
proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
proj <- addCombinedDims(proj, reducedDims = c("Harmony_ATAC", "Harmony_RNA"), name =  "Harmony_Combined")

#UMAPs
proj <- addUMAP(proj, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
#Add Clusters
proj <- addClusters(proj, reducedDims = "LSI_Combined", name = "Clusters", resolution = 0.2, force = TRUE)

#UMAPs
proj <- addUMAP(proj, reducedDims = "Harmony_ATAC", name = "UMAP_Harmony_ATAC", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "Harmony_RNA", name = "UMAP_Harmony_RNA", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "Harmony_Combined", name = "UMAP_Harmony_Combined", minDist = 0.8, force = TRUE)
#Add Clusters
proj <- addClusters(proj, reducedDims = "Harmony_Combined", name = "Harmony_Clusters", resolution = 0.2, force = TRUE)

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered_multiome_2')

#Plot Embedding
p1 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
#Save Plot
plotPDF(p1, p2, p3, name = "UMAP-scATAC-scRNA-Combined", addDOC = FALSE)

#Plot Embedding
p1 <- plotEmbedding(proj, name = "sample", embedding = "UMAP_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj, name = "sample", embedding = "UMAP_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = "sample", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
#Save Plot
plotPDF(p1, p2, p3, name = "UMAP-scATAC-scRNA-Combined-Batch", addDOC = FALSE)

#Plot Embedding
p1 <- plotEmbedding(proj, name = "condition", embedding = "UMAP_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj, name = "condition", embedding = "UMAP_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = "condition", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
#Save Plot
plotPDF(p1, p2, p3, name = "UMAP-scATAC-scRNA-Combined-Condition", addDOC = FALSE)

#Plot Embedding
p1 <- plotEmbedding(proj, name = "Harmony_Clusters", embedding = "UMAP_Harmony_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj, name = "Harmony_Clusters", embedding = "UMAP_Harmony_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = "Harmony_Clusters", embedding = "UMAP_Harmony_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
#Save Plot
plotPDF(p1, p2, p3, name = "UMAP-scATAC-scRNA-Combined-Harmony", addDOC = FALSE)

#Plot Embedding
p1 <- plotEmbedding(proj, name = "sample", embedding = "UMAP_Harmony_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj, name = "sample", embedding = "UMAP_Harmony_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = "sample", embedding = "UMAP_Harmony_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
#Save Plot
plotPDF(p1, p2, p3, name = "UMAP-scATAC-scRNA-Combined-Batch-Harmony", addDOC = FALSE)

#Plot Embedding
p1 <- plotEmbedding(proj, name = "celltype", embedding = "UMAP_Harmony_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
# p2 <- plotEmbedding(proj, name = "predicted.celltype.l1", embedding = "UMAP_Harmony_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
# p3 <- plotEmbedding(proj, name = "predicted.celltype.l2", embedding = "UMAP_Harmony_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
#Save Plot
plotPDF(p1, name = "UMAP-scATAC-scRNA-Combined-Predicted-Harmony", addDOC = FALSE)

data <- readRDS('../tcell_annotated_updated_2_conditions.RDS')
seuratUMAP <- function(name, seurat_name){
  proj@embeddings[[name]] <- proj@embeddings$UMAP_Combined
  table(rownames(proj@embeddings[[name]]$df) %in% gsub('_', '#', rownames(data@reductions[[seurat_name]]@cell.embeddings)))
  df <- data@reductions[[seurat_name]]@cell.embeddings
  rownames(df) <- gsub('_', '#', rownames(df))
  df <- df[rownames(proj@embeddings[[name]]$df),]
  colnames(df) <- colnames(proj@embeddings$UMAP_Combined$df)
  table(rownames(proj@embeddings[[name]]$df) == rownames(df))
  proj@embeddings[[name]]$df <- df
  return(proj)
}

proj <- seuratUMAP('wnn.umap', 'wnn.umap')
proj <- seuratUMAP('rna.umap', 'rna.umap')
proj <- seuratUMAP('adt.umap', 'adt.umap')
proj <- seuratUMAP('atac.umap', 'atac.umap')
proj <- seuratUMAP('wnn2.umap', 'wnn2.umap')
#Plot Embedding
p1 <- plotEmbedding(proj, name = "celltype", embedding = "wnn2.umap", size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj, name = "celltype", embedding = "wnn.umap", size = 1.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = "celltype", embedding = "rna.umap", size = 1.5, labelAsFactors=F, labelMeans=F)
p4 <- plotEmbedding(proj, name = "celltype", embedding = "adt.umap", size = 1.5, labelAsFactors=F, labelMeans=F)
p5 <- plotEmbedding(proj, name = "celltype", embedding = "atac.umap", size = 1.5, labelAsFactors=F, labelMeans=F)

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered_multiome_2')

#Save Plot
plotPDF(p1, p2, p3, p4, p5, name = "UMAP-Celltype", addDOC = FALSE)

#############
proj <- addImputeWeights(proj, reducedDims = 'Harmony_Combined')

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered_multiome_2')
