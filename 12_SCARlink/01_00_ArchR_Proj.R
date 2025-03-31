rm(list=ls())
gc()

.libPaths()
user_lib <- "/ihome/wchen/zhongli/R/x86_64-pc-linux-gnu-library/4.0"
.libPaths(c(user_lib, .libPaths()))

# Load Packages
library(data.table)
library(tidyverse)
library(dplyr)
library(Seurat)
library(ArchR)
library(Signac)
library(SeuratWrappers)
library(stringr)
library(parallel)

set.seed(1)

# Parameters
ArchR_PATH <- "/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Result/"

setwd(WORK_PATH)

addArchRThreads(threads = 15) 

paths <- c('Duerr_20210419_DOGMAseq_DIG_arc',
           'Duerr_20210610_DOGMAseq-1_arc', 'Duerr_20210610_DOGMAseq-2_arc',
           'Duerr_20210826_DOGMAseq-1_arc', 'Duerr_20210826_DOGMAseq-2_arc',
           'Duerr_20210831_DOGMAseq-1_arc', 'Duerr_20210831_DOGMAseq-2_arc',
           'Duerr_20220126_DOGMAseq_arc',
           'Duerr_20220215_DOGMAseq-1_arc', 'Duerr_20220215_DOGMAseq-2_arc',
           'Duerr_20220218_DOGMAseq-1_arc', 'Duerr_20220218_DOGMAseq-2_arc',
           'Duerr_20220524_DOGMAseq-1_arc', 'Duerr_20220524_DOGMAseq-2_arc',
           'Duerr_20220603_DOGMAseq-2_arc')

inputFiles <- paste0('/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/', paths, '/atac_fragments.tsv.gz')

names(inputFiles) <- c('04191', '06101', '06102', '08261', '08262',
                       '08311', '08312', '01261', '02151', '02152',
                       '02181', '02182', '05241', '05242', '06032')

addArchRGenome("hg38")

# system("cp /ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/DOGMA_filtered/ArrowFiles/* /ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/SCARlink/Result/ArrowFiles/")

# ArrowFiles <- createArrowFiles(inputFiles = inputFiles, sampleNames = names(inputFiles),
#                                minTSS = 0, #Dont set this too high because you can always increase later
#                                minFrags = 0, maxFrags = 1e+07,
#                                addTileMat = TRUE, addGeneScoreMat = TRUE)

ArrowFiles <- paste0(c('04191', '06101', '06102', '08261', '08262',
                       '08311', '08312', '01261', '02151', '02152',
                       '02181', '02182', '05241', '05242', '06032'), ".arrow")

proj <- ArchRProject(ArrowFiles = ArrowFiles,  outputDirectory = 'DOGMA_ATAC', copyArrows = FALSE)

getAvailableMatrices(proj)

meta <- read.csv('/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/tcell_annotated_updated.csv', row.names = 'X')

meta$cellNames <- rownames(meta)
meta$cellNames <- gsub('_', '#', meta$cellNames)

table(proj$cellNames %in% meta$cellNames)
proj <- proj[proj$cellNames %in% meta$cellNames,]
rownames(meta) <- meta$cellNames
meta <- meta[proj$cellNames,]
proj$condition <- meta$condition
proj$sample <- meta$sample
proj$celltype_updated <- meta$celltype_updated

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_ATAC')

# proj_list <- list()
# for (i in seq_along(inputFiles)) {
#   sample_name <- names(inputFiles)[i]
#   arrow_file <- paste0(sample_name, ".arrow")
#   
#   # Create ArchRProject for each Arrow file
#   proj <- ArchRProject(
#     ArrowFiles = arrow_file,
#     outputDirectory = paste0("ArchRProject_", sample_name),
#     copyArrows = FALSE)
#   
#   # Add the TileMatrix
#   proj <- addTileMatrix(input = proj, binarize = FALSE, tileSize = 500, force = TRUE)
#   
#   tm <- getMatrixFromProject(ArchRProj = proj,
#                              useMatrix = 'TileMatrix', binarize = FALSE)
#   
#   proj <- saveArchRProject(ArchRProj = proj, outputDirectory = paste0("SavedProject_", sample_name))
#   proj_list[[sample_name]] <- proj
#   print(i)
# }

