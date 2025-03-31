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
library(parallel)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')

load('../../Duerr_20210419_DOGMAseq_DIG_arc/Analysis/data.RData')
data1 <- data
rm(data)
load('../../Duerr_20210610_DOGMAseq-1_arc/Analysis/data.RData')
data2 <- data
rm(data)
load('../../Duerr_20210610_DOGMAseq-2_arc/Analysis/data.RData')
data3 <- data
rm(data)
load('../../Duerr_20210826_DOGMAseq-1_arc/Analysis/data.RData')
data4 <- data
rm(data)
load('../../Duerr_20210826_DOGMAseq-2_arc/Analysis/data.RData')
data5 <- data
rm(data)
load('../../Duerr_20210831_DOGMAseq-1_arc/Analysis/data.RData')
data6 <- data
rm(data)
load('../../Duerr_20210831_DOGMAseq-2_arc/Analysis/data.RData')
data7 <- data
rm(data)
load('../../Duerr_20220126_DOGMAseq_arc/Analysis/data.RData')
data8 <- data
rm(data)
load('../../Duerr_20220215_DOGMAseq-1_arc/Analysis/data.RData')
data9 <- data
rm(data)
load('../../Duerr_20220215_DOGMAseq-2_arc/Analysis/data.RData')
data10 <- data
rm(data)
load('../../Duerr_20220218_DOGMAseq-1_arc/Analysis/data.RData')
data11 <- data
rm(data)
load('../../Duerr_20220218_DOGMAseq-2_arc/Analysis/data.RData')
data12 <- data
rm(data)
load('../../Duerr_20220524_DOGMAseq-1_arc/Analysis/data.RData')
data13 <- data
rm(data)
load('../../Duerr_20220524_DOGMAseq-2_arc/Analysis/data.RData')
data14 <- data
rm(data)
load('../../Duerr_20220603_DOGMAseq-2_arc/Analysis/data.RData')
data15 <- data
rm(data)

data_idx <- paste0('data', 1:15)
id_idx <- c('04191',
            '06101', '06102',
            '08261', '08262',
            '08311', '08312',
            '01261',
            '02151', '02152',
            '02181', '02182',
            '05241', '05242',
            '06032')

# for(i in 1:15){
#   assign(data_idx[i], RenameCells(get(data_idx[i]), add.cell.id = id_idx[i]))
# }

extract_matrix <- function(data, type, id){
  if(type == 'rna'){
    data <- data@assays$RNA@counts
    colnames(data) <- paste(id, colnames(data), sep = '_')
  }
  if(type == 'adt'){
    data <- data@assays$ADT@counts
    colnames(data) <- paste(id, colnames(data), sep = '_')
  }
  if(type == 'meta'){
    data <- data@meta.data
    rownames(data) <- paste(id, rownames(data), sep = '_')
  }
  return(data)
}

rna_matx <- mclapply(1:15, function(x){extract_matrix(get(data_idx[x]), 'rna', id_idx[x])}, mc.cores = 15)
rna <- Reduce(cbind, rna_matx)

adt_matx <- mclapply(1:15, function(x){extract_matrix(get(data_idx[x]), 'adt', id_idx[x])}, mc.cores = 15)
adt <- Reduce(cbind, adt_matx)

meta_matx <- mclapply(1:15, function(x){extract_matrix(get(data_idx[x]), 'meta', id_idx[x])}, mc.cores = 15)
meta <- Reduce(rbind, meta_matx)

extract_fragments <- function(h5_path, data, id){
  fragments <- CreateFragmentObject(
    path = paste0(h5_path,'/atac_fragments.tsv.gz'),
    cells = colnames(data),
    validate.fragments = FALSE
  )
  names(fragments@cells) <- paste(id, names(fragments@cells), sep = '_')
  return(fragments)
}

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

fragments <- mclapply(1:15, function(x){extract_fragments(paste0('../../', paths[x]), get(data_idx[x]), id_idx[x])}, mc.cores = 15)

data <- CreateSeuratObject(
  counts = rna,
  assay = "RNA",
  meta.data = meta
)

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data[['ADT']] <- CreateAssayObject(adt)

save(data, fragments, file = '../data/data_fragments.RData')
barcodes <- as.character(colnames(data))
save(barcodes, file = '../data/barcodes.RData')

load('../data/data_fragments.RData')

atac <- matrix(1, nrow = 2, ncol = dim(data)[2])
rownames(atac) <- c('1-1-2', '1-3-4')
colnames(atac) <- colnames(data)
data[["peaks"]] <- CreateChromatinAssay(
  counts = atac,
  fragments = fragments,
  validate.fragments = F
)

DefaultAssay(data) <- "peaks"

plan('multisession')
options(future.globals.maxSize= 36*1024^3)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

peaks <- CallPeaks(data, macs2.path = "/ihome/wchen/zhongli/.conda/envs/MACS2/bin/macs2")
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

save(peaks, file = '../data/peaks.RData')

load('../data/peaks.RData')

macs2_counts <- FeatureMatrix(
  fragments = fragments,
  features = peaks,
  cells = colnames(data)
)

save(macs2_counts, file = '../data/macs2_counts.RData')

# create a new assay using the MACS2 peak set and add it to the Seurat object
data[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragments,
  annotation = annotations
)

save(data, file = '../output/data_macs2.RData')