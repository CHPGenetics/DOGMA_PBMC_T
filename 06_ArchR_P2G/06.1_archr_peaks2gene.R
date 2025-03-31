library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
set.seed(1)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/')

addArchRThreads(threads = 64)

proj <- loadArchRProject('DOGMA_filtered_multiome/')

proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "Harmony_Combined",
  useMatrix = "GeneExpressionMatrix",
  maxDist = 1000000
)

saveArchRProject(proj, 'DOGMA_filtered_multiome')

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.4,
  resolution = 1,
  returnLoops = FALSE
)

saveRDS(p2g, file = 'DOGMA_filtered_multiome/p2g.RDS')

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0,
  FDRCutOff= 0.05,
  resolution = 1,
  returnLoops = FALSE
)

saveRDS(p2g, file = 'DOGMA_filtered_multiome/p2g_fdr_0.05.RDS')

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0,
  FDRCutOff= 0.05,
  resolution = 1
)

saveRDS(p2g, file = 'DOGMA_filtered_multiome/p2g_fdr_0.05_loop.RDS')

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0,
  FDRCutOff= 0.01,
  resolution = 1,
  returnLoops = FALSE
)

saveRDS(p2g, file = 'DOGMA_filtered_multiome/p2g_fdr_0.01.RDS')

p <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "celltype_id", corCutOff = 0,
                          FDRCutOff= 0.05)
plotPDF(p, name = "Peak2Gene-Heatmap", width = 16, height = 16, ArchRProj = proj, addDOC = FALSE)

p <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "celltype_id", corCutOff = 0,
                          FDRCutOff= 0.05, returnMatrices = T)
save(p, file = 'DOGMA_filtered_multiome/p2g_heatmap_matrix.RData')