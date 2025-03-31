library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
set.seed(1)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/')

addArchRThreads(threads = 64)

proj <- loadArchRProject('DOGMA_filtered_multiome_2/')

p2g_analysis <- function(celltype, condition){
  proj_sub <- proj[proj$celltype_id == celltype,]
  proj_sub <- proj_sub[proj_sub$condition == condition,]
  
  dir.create(file.path('DOGMA_filtered_multiome_2/', celltype))
  dir.create(file.path(paste0('DOGMA_filtered_multiome_2/', celltype), condition))
  
  proj_sub <- addPeak2GeneLinks(
    ArchRProj = proj_sub,
    reducedDims = "Harmony_Combined",
    useMatrix = "GeneExpressionMatrix",
    maxDist = 1000000
  )
  
  p2g <- getPeak2GeneLinks(
    ArchRProj = proj_sub,
    corCutOff = -1,
    FDRCutOff= 0.001,
    resolution = 1,
    returnLoops = FALSE
  )

  saveRDS(p2g, file = paste0('DOGMA_filtered_multiome_2/', celltype, '/', condition, '/p2g_fdr_0.001.RDS'))
  
  p2g <- getPeak2GeneLinks(
    ArchRProj = proj_sub,
    corCutOff = -1,
    FDRCutOff= 0.001,
    resolution = 1
  )
  
  saveRDS(p2g, file = paste0('DOGMA_filtered_multiome_2/', celltype, '/', condition, '/p2g_fdr_0.001_loop.RDS'))
}

for(i in names(table(proj$celltype_id))){
    for(j in names(table(proj$condition))){
        p2g_analysis(i,j)
    }
}