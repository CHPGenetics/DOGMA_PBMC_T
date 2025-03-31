library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
set.seed(1)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/')

addArchRThreads(threads = 64)

proj <- loadArchRProject('DOGMA_filtered/')

paths <- c('Duerr_20210419_DOGMAseq_DIG_arc',
           'Duerr_20210610_DOGMAseq-1_arc', 'Duerr_20210610_DOGMAseq-2_arc',
           'Duerr_20210826_DOGMAseq-1_arc', 'Duerr_20210826_DOGMAseq-2_arc',
           'Duerr_20210831_DOGMAseq-1_arc', 'Duerr_20210831_DOGMAseq-2_arc',
           'Duerr_20220126_DOGMAseq_arc',
           'Duerr_20220215_DOGMAseq-1_arc', 'Duerr_20220215_DOGMAseq-2_arc',
           'Duerr_20220218_DOGMAseq-1_arc', 'Duerr_20220218_DOGMAseq-2_arc',
           'Duerr_20220524_DOGMAseq-1_arc', 'Duerr_20220524_DOGMAseq-2_arc',
           'Duerr_20220603_DOGMAseq-2_arc'
)

inputFiles <- paste0('../../../', paths, '/filtered_feature_bc_matrix.h5')

seRNA <- import10xFeatureMatrix(
  input = inputFiles,
  names = c('04191',
            '06101', '06102',
            '08261', '08262',
            '08311', '08312',
            '01261',
            '02151', '02152',
            '02181', '02182',
            '05241', '05242',
            '06032')
)
rna <- Reduce(cbind, seRNA)
proj <- addGeneExpressionMatrix(input = proj, seRNA = rna, force = TRUE)

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered_multiome')

# rna <- rna[,proj$cellNames]
# load('../../output/data_harmony.RData')
# adt <- as.SingleCellExperiment(data, assay = 'ADT')
# colnames(adt) <- gsub('0419_', 'Duerr_20210419_DOGMAseq#', colnames(adt))
# colnames(adt) <- gsub('0610_', 'Duerr_20210610_DOGMAseq#', colnames(adt))
# colnames(adt) <- gsub('0826_', 'Duerr_20210826_DOGMAseq#', colnames(adt))
# colnames(adt) <- gsub('0831_', 'Duerr_20210831_DOGMAseq#', colnames(adt))
# adt <- adt[,colnames(rna)]
# 
# antibody <- read.xlsx('TotalSeq_A_Human_Universal_Cocktail_v1_163_Antibodies_399907_Barcodes.xlsx')
# 
# adt <- assay(adt, 'counts')
# adt <- adt[!is.na(antibody$Gene.Name),]
# antibody <- antibody[!is.na(antibody$Gene.Name),]
# rownames(adt) <- antibody$Gene.Name
# 
# adt.se <- rna[rownames(adt),]
# assay(adt.se) <- adt
# adt <- adt.se
# 
# proj <- addFeatureMatrix(input = proj, features = rowRanges(adt), force = TRUE, matrixName = "ProteinMatrix")

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
  name = "LSI_ATAC"
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
  name = "LSI_RNA"
)

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "LSI_ATAC",
  name = "Harmony_ATAC",
  groupBy = c("sample")
)

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "LSI_RNA",
  name = "Harmony_RNA",
  groupBy = c("sample")
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

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered_multiome')

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



