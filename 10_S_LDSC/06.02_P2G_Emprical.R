rm(list=ls())
gc()

# Load Packages
library(Seurat)
library(Signac)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ArchR)
library(parallel)
library(readxl)

# Parameters
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/04_Global_Analysis/04_P2G_Shuffle/"

# Function
sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
  // [[Rcpp::export]]
  Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
    
    if(X.ncol() != Y.ncol()){
      stop("Columns of Matrix X and Y must be equal length!");
    }

    if(max(idxX) > X.nrow()){
      stop("Idx X greater than nrow of Matrix X");
    }

    if(max(idxY) > Y.nrow()){
      stop("Idx Y greater than nrow of Matrix Y");
    }

    // Transpose Matrices
    X = transpose(X);
    Y = transpose(Y);
    
    const int nx = X.ncol();
    const int ny = Y.ncol();

    // Centering the matrices
    for (int j = 0; j < nx; ++j) {
      X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
    }

    for (int j = 0; j < ny; ++j) {
      Y(Rcpp::_, j) = Y(Rcpp::_, j) - Rcpp::mean(Y(Rcpp::_, j));
    }

    // Compute 1 over the sample standard deviation
    Rcpp::NumericVector inv_sqrt_ss_X(nx);
    for (int i = 0; i < nx; ++i) {
      inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum( X(Rcpp::_, i) * X(Rcpp::_, i) ));
    }

    Rcpp::NumericVector inv_sqrt_ss_Y(ny);
    for (int i = 0; i < ny; ++i) {
      inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum( Y(Rcpp::_, i) * Y(Rcpp::_, i) ));
    }

    //Calculate Correlations
    const int n = idxX.size();
    Rcpp::NumericVector cor(n);
    for(int k = 0; k < n; k++){
      cor[k] = Rcpp::sum( X(Rcpp::_, idxX[k] - 1) * Y(Rcpp::_, idxY[k] - 1) ) * inv_sqrt_ss_X(idxX[k] - 1) * inv_sqrt_ss_Y(idxY[k] - 1);
    } 

    return(cor);

  }'
)

.getNullCorrelations <- function(seA, seB, o, n){
  
  o$seq <- as.character(seqnames(seA)[o$peak_idx])
  
  nullCor <- lapply(seq_along(unique(o$seq)), function(i){
    
    #Get chr from olist
    chri <- unique(o$seq)[i]
    message(chri, " ", appendLF = FALSE)
    
    #Randomly get n seA
    id <- which(as.character(seqnames(seA)) != chri)
    if(length(id) > n){
      transAidx <- sample(id, n)
    }else{
      transAidx <- id
    }
    
    #Calculate Correlations
    grid <- expand.grid(transAidx, unique(o[o$seq==chri,]$gene_idx))
    
    idxA <- unique(grid[,1])
    idxB <- unique(grid[,2])
    
    seSubA <- seA[idxA]
    seSubB <- seB[idxB]
    
    grid[,3] <- match(grid[,1], idxA)
    grid[,4] <- match(grid[,2], idxB)
    
    colnames(grid) <- c("A", "B")
    out <- rowCorCpp(grid[,3], grid[,4], as.matrix(assay(seSubA)), as.matrix(assay(seSubB)))
    out <- na.omit(out)
    
    return(out)
    
  }) %>% SimpleList
  message("")
  
  summaryDF <- lapply(nullCor, function(x){
    data.frame(mean = mean(x), sd = sd(x), median = median(x), n = length(x))
  }) %>% Reduce("rbind",.)
  
  return(list(summaryDF, unlist(nullCor)))
  
}

.getQuantiles <- function(v = NULL, len = length(v)){
  if(length(v) < len){
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  }else{
    v2 <- v
  }
  p <- trunc(rank(v2))/length(v2)
  if(length(v) < len){
    p <- p[seq_along(v)]
  }
  return(p)
}

peak2gene_emprical <- function(peak.mat = NULL,
                               rna.mat = NULL,
                               max.dist = 250000){
  
  # Gene Annotation
  gene_anno <- geneAnnoHg38
  genes <- gene_anno$genes
  gene.use <- intersect(elementMetadata(genes)[, "symbol"], rownames(rna.mat))
  genes <- genes[elementMetadata(genes)[, "symbol"] %in% gene.use]
  rna.mat <- rna.mat[gene.use, ]
  
  gene_start <- ifelse(genes@strand == "+",
                       genes@ranges@start,
                       genes@ranges@start + genes@ranges@width - 1)
  
  genes <- GRanges(genes@seqnames,
                   ranges = IRanges(gene_start, width = 1),
                   name = genes$symbol,
                   gene_id = genes$gene_id,
                   strand = genes@strand)
  
  seRNA <- SummarizedExperiment(assays = SimpleList(RNA = rna.mat), rowRanges = genes)
  
  # ATAC-seq data
  df_peak <- stringr::str_split_fixed(rownames(peak.mat), "-", 3)
  
  peakSet <- GRanges(df_peak[, 1],
                     IRanges(start = as.numeric(df_peak[, 2]),
                             end = as.numeric(df_peak[, 3])))
  
  seATAC <- SummarizedExperiment(assays = SimpleList(ATAC = peak.mat), rowRanges = peakSet)
  
  # Putative Peak-to-gene
  o <- data.frame(findOverlaps(
    resize(seRNA, 2 * max.dist + 1, "center"),
    resize(rowRanges(seATAC), 1, "center"),
    ignore.strand = TRUE))
  o$distance <- IRanges::distance(rowRanges(seRNA)[o[, 1]],
                                  rowRanges(seATAC)[o[, 2]])
  colnames(o) <- c("gene_idx", "peak_idx", "distance")
  
  ## Null Correlations
  nullCor <- .getNullCorrelations(seATAC, seRNA, o, 1000)
  
  rm(.Random.seed, envir=globalenv())
  n <- floor(runif(1, min=0, max=1001))
  set.seed(n)
  seATAC.mat <- assay(seATAC)
  seATAC.mat <- seATAC.mat[,sample(ncol(seATAC.mat))] # Permute ATAC aggregates
  
  o$Correlation <- rowCorCpp(as.integer(o$peak_idx), as.integer(o$gene_idx), 
                             as.matrix(seATAC.mat), as.matrix(assay(seRNA)))
  
  df <- rowRanges(seATAC)[o$peak_idx, ]
  
  o$gene <- rowData(seRNA)[o$gene_idx, ]$name
  o$peak <- paste0(df@seqnames, "-", as.data.frame(df@ranges)$start, 
                   "-", as.data.frame(df@ranges)$end)
  
  # Correlation
  df <- rowRanges(seATAC)[o$peak_idx, ]
  
  o$gene <- rowData(seRNA)[o$gene_idx, ]$name
  o$peak <- paste0(df@seqnames, "-", as.data.frame(df@ranges)$start, "-",
                   as.data.frame(df@ranges)$end)
  
  ## compute correlation
  o$Correlation <- rowCorCpp(as.integer(o$peak_idx),
                             as.integer(o$gene_idx),
                             as.matrix(assay(seATAC)),
                             as.matrix(assay(seRNA)))
  
  ## compute p-value
  o$VarAssayA <- .getQuantiles(matrixStats::rowVars(as.matrix(seATAC.mat)))[o$peak_idx]
  o$VarAssayB <- .getQuantiles(matrixStats::rowVars(as.matrix(assay(seRNA))))[o$gene_idx]
  o$TStat <- (o$Correlation / sqrt((1-o$Correlation^2)/(ncol(seATAC.mat)-2))) #T-statistic P-value
  o$Pval <- 2*pt(-abs(o$TStat), ncol(seATAC.mat) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  o$EmpPval <- 2*pnorm(-abs(((o$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))
  o$EmpFDR <- p.adjust(o$EmpPval, method = "fdr")
  
  return(o)
}

# Load data
rna.mat.raw <- readRDS(paste0(WORK_PATH, "rna.rds"))
peak.mat.raw <- readRDS(paste0(WORK_PATH, "peak.rds"))
meta <- readRDS(paste0(WORK_PATH, "meta.rds"))

# P2G in specific cell type
Result <- data.frame()
for(celltype in unique(meta$celltype_updated)){
  
  barcode <- rownames(meta[which(meta$celltype_updated == celltype),])
  rna.mat <- rna.mat.raw[,barcode]
  peak.mat <- peak.mat.raw[,barcode]
  
  p2g_result <- peak2gene_emprical(peak.mat = peak.mat,
                                   rna.mat = rna.mat,
                                   max.dist = 250000)
  p2g_result$cell_type <- celltype
  saveRDS(p2g_result, paste0(WORK_PATH, "P2G/Emprical/P2G_", celltype, ".rds"))
  
  Result <- rbind(Result, p2g_result)
  print(celltype)
}
saveRDS(Result, paste0(WORK_PATH, "P2G/Emprical/All_P2G.rds"))
