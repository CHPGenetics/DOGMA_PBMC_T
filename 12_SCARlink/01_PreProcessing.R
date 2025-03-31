rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(dplyr)
library(Seurat)
library(ArchR)
library(Signac)
library(SeuratWrappers)

addArchRThreads(threads = 64) 

# Parameters
ArchR_PATH <- "/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/"
WORK_PATH <- "/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Result/"

setwd(WORK_PATH)

# tcell_dogma <- readRDS("/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/tcell_annotated_updated.RDS")
# 
# scrna.object <- CreateSeuratObject(counts = tcell_dogma[["RNA"]]@counts, 
#                                    meta.data = tcell_dogma@meta.data)
# scrna.object <- AddMetaData(scrna.object, metadata = tcell_dogma@meta.data)
# 
# scatac.object <- loadArchRProject(paste0(WORK_PATH, 'DOGMA_ATAC/')) # ArchR project
# 
# head(colnames(scrna.object))
# head(scatac.object$cellNames)
# scrna.object <- RenameCells(scrna.object, new.names=gsub("_", "#", colnames(scrna.object)))
# saveRDS(scrna.object, paste0(WORK_PATH, "dogma_scrna.rds"))

scrna.object <- readRDS(paste0(WORK_PATH, "dogma_scrna.rds"))
scatac.object <- loadArchRProject(paste0(WORK_PATH, 'DOGMA_ATAC/'))

rna_cells <- colnames(scrna.object)
atac_cells <- scatac.object$cellNames

set.seed(1)
selected_cells <- sample(rna_cells, size = 100000, replace = FALSE)

scatac.object <- subsetArchRProject(ArchRProj = scatac.object, cells = selected_cells, 
                                    outputDirectory = "DOGMA_ATAC_Subset")
scrna.object <- subset(scrna.object, cells = selected_cells)

# Order cells
rna_cells <- colnames(scrna.object)
atac_cells <- scatac.object$cellNames
idxMatch <- match(rna_cells, atac_cells)
scatac.object <- scatac.object[idxMatch, ]

# Relative counts normalization scRNA-seq count matrix
scrna.object <- FindVariableFeatures(scrna.object, selection.method = "vst", nfeatures = 5000)
scrna.object <- NormalizeData(scrna.object, normalization.method = "RC")
scrna.object$celltype_updated <- as.character(scrna.object$celltype_updated)
saveRDS(scrna.object, paste0(WORK_PATH, "dogma_scrna_subset.rds"))

addTileMatrix(input = scatac.object, binarize = FALSE, tileSize = 500, force = TRUE)

saveArchRProject(scatac.object, paste0(WORK_PATH, "DOGMA_ATAC_Subset"))

tm <- getMatrixFromProject(ArchRProj = scatac.object,
                           useMatrix = 'TileMatrix', binarize = FALSE)

scatac.object <- addIterativeLSI(ArchRProj = scatac.object,
                                 useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2, 
                                 clusterParams = list(resolution = c(0.2), 
                                                      sampleCells = 10000, 
                                                      n.start = 10), 
                                 varFeatures = 25000, dimsToUse = 1:30)

saveArchRProject(scatac.object, paste0(WORK_PATH, "DOGMA_ATAC_Subset"))

# Split into four conditions
scrna.object <- readRDS(paste0(WORK_PATH, "dogma_scrna_subset.rds"))
scatac.object <- loadArchRProject(paste0(WORK_PATH, 'DOGMA_ATAC_Subset/'))
## scRNA
scrna.object <- readRDS(paste0(WORK_PATH, "dogma_scrna_subset.rds"))
scrna.object <- FindVariableFeatures(scrna.object, selection.method = "vst", nfeatures = 5000)
scrna.object <- NormalizeData(scrna.object, normalization.method = "RC")
scrna.object$celltype_updated <- as.character(scrna.object$celltype_updated)

for(i in c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2", 
           "Act_IL1B_IL23_TGFB", "Act_IL1B_IL23_PGE2_TGFB")){
  
  scrna.object_sub <- subset(scrna.object, subset = condition == i)
  
  # Select the cells
  set.seed(2024)
  selected_cells <- sample(colnames(scrna.object_sub), size = 10000, replace = FALSE)
  scatac.object_sub <- subsetArchRProject(ArchRProj = scatac.object, cells = selected_cells, 
                                          outputDirectory = paste0(i, "/DOGMA_ATAC_Subset"))
  scrna.object_sub <- subset(scrna.object_sub, cells = selected_cells)
  
  # Order
  rna_cells <- colnames(scrna.object_sub)
  atac_cells <- scatac.object_sub$cellNames
  idxMatch <- match(rna_cells, atac_cells)
  scatac.object_sub <- scatac.object_sub[idxMatch, ]
  
  # Save
  saveRDS(scrna.object_sub, paste0(WORK_PATH, i, "/dogma_scrna_subset.rds"))
  saveArchRProject(scatac.object_sub, paste0(WORK_PATH, i, "/DOGMA_ATAC_Subset"))
  
  print(i)
}

# All
# Select the cells
set.seed(2024)
selected_cells <- sample(colnames(scrna.object), size = 10000, replace = FALSE)
scatac.object_sub <- subsetArchRProject(ArchRProj = scatac.object, cells = selected_cells, 
                                        outputDirectory = paste0("All/DOGMA_ATAC_Subset"))
scrna.object_sub <- subset(scrna.object, cells = selected_cells)

# Order
rna_cells <- colnames(scrna.object_sub)
atac_cells <- scatac.object_sub$cellNames
idxMatch <- match(rna_cells, atac_cells)
scatac.object_sub <- scatac.object_sub[idxMatch, ]

# Save
saveRDS(scrna.object_sub, paste0(WORK_PATH, "All/dogma_scrna_subset.rds"))
saveArchRProject(scatac.object_sub, paste0(WORK_PATH, "All/DOGMA_ATAC_Subset"))

# 50k
set.seed(2024)
selected_cells <- sample(colnames(scrna.object), size = 50000, replace = FALSE)
scatac.object_sub <- subsetArchRProject(ArchRProj = scatac.object, cells = selected_cells, 
                                        outputDirectory = paste0("New_50k/DOGMA_ATAC_Subset"))
scrna.object_sub <- subset(scrna.object, cells = selected_cells)

# Order
rna_cells <- colnames(scrna.object_sub)
atac_cells <- scatac.object_sub$cellNames
idxMatch <- match(rna_cells, atac_cells)
scatac.object_sub <- scatac.object_sub[idxMatch, ]

# Save
saveRDS(scrna.object_sub, paste0(WORK_PATH, "New_50k/dogma_scrna_subset.rds"))
saveArchRProject(scatac.object_sub, paste0(WORK_PATH, "New_50k/DOGMA_ATAC_Subset"))



