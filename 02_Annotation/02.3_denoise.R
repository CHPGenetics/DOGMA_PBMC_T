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
library(dsb)
library(celda)
library(parallel)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')

load('../output/data_harmony_qc_tcell_updated.RData')
barcodes <- as.character(colnames(data))
save(barcodes, file = '../data/tcell_barcodes.RData')
rm(data)

load('../data/tcell_barcodes.RData')
df <- as.data.frame(str_split_fixed(barcodes, '_', 2))
df$sample <- df$V1
samples <- c('Duerr_20210419_DOGMAseq_DIG',
           'Duerr_20210610_DOGMAseq-1', 'Duerr_20210610_DOGMAseq-2',
           'Duerr_20210826_DOGMAseq-1', 'Duerr_20210826_DOGMAseq-2',
           'Duerr_20210831_DOGMAseq-1', 'Duerr_20210831_DOGMAseq-2',
           'Duerr_20220126_DOGMAseq',
           'Duerr_20220215_DOGMAseq-1', 'Duerr_20220215_DOGMAseq-2',
           'Duerr_20220218_DOGMAseq-1', 'Duerr_20220218_DOGMAseq-2',
           'Duerr_20220524_DOGMAseq-1', 'Duerr_20220524_DOGMAseq-2',
           'Duerr_20220603_DOGMAseq-2'
)
index <- c('04191',
            '06101', '06102',
            '08261', '08262',
            '08311', '08312',
            '01261',
            '02151', '02152',
            '02181', '02182',
            '05241', '05242',
            '06032')

for (i in 1:length(index)){
  df$sample <- gsub(index[i], samples[i], df$sample)
}

denoise_adt <- function(dataset){
  h5_path <- paste0(dataset, '_arc')
  metadata <- read.csv(
    file = paste0('../../', h5_path,'/per_barcode_metrics.csv'),
    header = TRUE,
    row.names = 1
  )
  metadata <- metadata[metadata$is_cell == 1,]
  
  h5_path2 <- dataset
  matx <- readMM(paste0('../../ADT/', h5_path2, '/featurecounts/featurecounts.mtx'))
  rownames(matx) <- fread(paste0('../../ADT/', h5_path2, '/featurecounts/featurecounts.barcodes.txt'), header = FALSE)[[1]]
  colnames(matx) <- fread(paste0('../../ADT/', h5_path2, '/featurecounts/featurecounts.genes.txt'), header = FALSE)[[1]]
  
  matx <- t(matx)
  colnames(matx) <- paste(colnames(matx), 1, sep = '-')
  rownames(matx) <- str_split_fixed(str_split_fixed(rownames(matx), '-A0', 2)[,1], '-A1', 2)[,1]
  
  data <- df[df$sample == dataset,]
  table(data$V2 %in% colnames(matx))                       
  # 
  # TRUE 
  # 20000 
  cell.adt.raw <- matx[,data$V2]
  backgroud <- setdiff(colnames(matx), rownames(metadata))
  background.adt.mtx <- matx[,backgroud]
  
  isotype.controls <- rownames(cell.adt.raw)[grep('Ctrl', rownames(cell.adt.raw))]
  
  cells.dsb.norm <- DSBNormalizeProtein(
    cell_protein_matrix = cell.adt.raw, 
    empty_drop_matrix = background.adt.mtx, 
    denoise.counts = TRUE, 
    use.isotype.control = TRUE, 
    isotype.control.name.vec = isotype.controls, 
    return.stats = TRUE
  )
  
  # counts <- Read10X_h5(paste0('../../', h5_path,'/raw_feature_bc_matrix.h5'))
  # rna_counts <- counts$`Gene Expression`
  # rm(counts)
  # 
  # rna_counts <- CreateSeuratObject(rna_counts)
  # rna_counts <- as.SingleCellExperiment(rna_counts)
  # 
  # cell.rna.raw <- rna_counts[,data$V2]
  # backgroud <- setdiff(colnames(rna_counts), rownames(metadata))
  # background.rna.mtx <- rna_counts[,backgroud]
  # 
  # sce <- decontX(cell.rna.raw, background = background.rna.mtx)
  # decontx <- decontXcounts(sce)
  saveRDS(cells.dsb.norm, file = paste0('../output/denoise/ADT/', dataset, '.RDS'))
  
  matx <- cells.dsb.norm$dsb_normalized_matrix
  colnames(matx) <- paste(names(table(data$V1)), colnames(matx), sep = '_')
  return(matx)
}

denoise_rna <- function(dataset){
  h5_path <- paste0(dataset, '_arc')
  metadata <- read.csv(
    file = paste0('../../', h5_path,'/per_barcode_metrics.csv'),
    header = TRUE,
    row.names = 1
  )
  metadata <- metadata[metadata$is_cell == 1,]
  
  # h5_path2 <- dataset
  # matx <- readMM(paste0('../../ADT/', h5_path2, '/featurecounts/featurecounts.mtx'))
  # rownames(matx) <- fread(paste0('../../ADT/', h5_path2, '/featurecounts/featurecounts.barcodes.txt'), header = FALSE)[[1]]
  # colnames(matx) <- fread(paste0('../../ADT/', h5_path2, '/featurecounts/featurecounts.genes.txt'), header = FALSE)[[1]]
  # 
  # matx <- t(matx)
  # colnames(matx) <- paste(colnames(matx), 1, sep = '-')
  # rownames(matx) <- str_split_fixed(str_split_fixed(rownames(matx), '-A0', 2)[,1], '-A1', 2)[,1]
  
  data <- df[df$sample == dataset,]
  # table(data$V2 %in% colnames(matx))                       
  # # 
  # # TRUE 
  # # 20000 
  # cell.adt.raw <- matx[,data$V2]
  # backgroud <- setdiff(colnames(matx), rownames(metadata))
  # background.adt.mtx <- matx[,backgroud]
  # 
  # isotype.controls <- rownames(cell.adt.raw)[grep('Ctrl', rownames(cell.adt.raw))]
  # 
  # cells.dsb.norm <- DSBNormalizeProtein(
  #   cell_protein_matrix = cell.adt.raw, 
  #   empty_drop_matrix = background.adt.mtx, 
  #   denoise.counts = TRUE, 
  #   use.isotype.control = TRUE, 
  #   isotype.control.name.vec = isotype.controls, 
  #   return.stats = TRUE
  # )
  
  counts <- Read10X_h5(paste0('../../', h5_path,'/raw_feature_bc_matrix.h5'))
  rna_counts <- counts$`Gene Expression`
  rm(counts)
  
  rna_counts <- CreateSeuratObject(rna_counts)
  rna_counts <- as.SingleCellExperiment(rna_counts)
  
  cell.rna.raw <- rna_counts[,data$V2]
  backgroud <- setdiff(colnames(rna_counts), rownames(metadata))
  background.rna.mtx <- rna_counts[,backgroud]
  
  sce <- decontX(cell.rna.raw, background = background.rna.mtx)
  decontx <- decontXcounts(sce)
  
  saveRDS(sce, file = paste0('../output/denoise/RNA/', dataset, '.RDS'))
  
  matx <- decontx
  colnames(matx) <- paste(names(table(data$V1)), colnames(matx), sep = '_')
  return(matx)
}

adt_matx <- mclapply(samples, function(x){denoise_adt(x)}, mc.cores = 15)
adt <- Reduce(cbind, adt_matx)

rna_matx <- mclapply(samples, function(x){denoise_rna(x)}, mc.cores = 15)
rna <- Reduce(cbind, rna_matx)

save(rna, adt, file = '../output/denoise/RNA_ADT.RData')