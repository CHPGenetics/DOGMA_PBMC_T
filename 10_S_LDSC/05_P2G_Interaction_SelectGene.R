rm(list=ls())
gc()

# Load Packages
library(Seurat)
library(Signac)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ArchR)
library(lme4)
library(lmerTest)
library(parallel)
library(readxl)

# Parameters
DATA_PATH <- "/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/04_Global_Analysis/03_P2G/"

# Function
# The function may have some problems
peak2gene_glmm <- function(peak.mat = NULL,
                           rna.mat = NULL,
                           max.dist = 250000,
                           condition.mat = NULL){
  
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
                   ranges = IRanges(gene_start,width = 1),
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
  
  df <- rowRanges(seATAC)[o$peak_idx, ]
  
  o$gene <- rowData(seRNA)[o$gene_idx, ]$name
  o$peak <- paste0(df@seqnames, "-", as.data.frame(df@ranges)$start, 
                   "-", as.data.frame(df@ranges)$end)
  
  # Mixed model
  results <- lapply(1:nrow(o), function(i) {
    gene_idx <- o$gene_idx[i]
    peak_idx <- o$peak_idx[i]
    
    geneValue <- assay(seRNA)[gene_idx, ] %>% as.numeric()
    peakValue <- assay(seATAC)[peak_idx, ] %>% as.numeric()
    
    p2g_data <- data.frame(gene = geneValue,
                           peak = peakValue,
                           condition = condition.mat$Condition)
    p2g_data <- p2g_data %>%
      filter(!(gene == 0 & peak == 0))
    
    if(nrow(p2g_data) > 30 & sum(p2g_data$gene) > 0  & sum(p2g_data$peak) > 0){
      
      model_glmm <- tryCatch({
        glmer(gene ~ peak + (peak | condition), 
              data = p2g_data, family = poisson())
      }, error = function(e) {
        
        if(grepl("pwrssUpdate did not converge", e$message)) {
          
          return(NULL)
        } else {
          
          stop(e)
        }
      })
      
      if(!is.null(model_glmm)){
        
        summary_model_glmm <- summary(model_glmm)
        
        ## Wald
        coef_est <- summary_model_glmm$coefficients["peak","Estimate"]
        std_error <- summary_model_glmm$coefficients["peak","Std. Error"]
        z_value <- coef_est / std_error
        p_wald <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
        
        ## LRTs
        full_model <- lmer(gene ~ peak + (peak | condition), data = p2g_data)
        reduced_model <- lmer(gene ~ 1 + (peak | condition), data = p2g_data)
        lrt <- anova(reduced_model, full_model)
        p_lrt <- lrt["full_model", "Pr(>Chisq)"]
        
      } else{
        
        p_wald = NA
        stat_wald = NA
        p_lrt = NA
        coef_est = NA
      }
      
      # GLM
      model_glm <- glm(gene ~ peak * condition, 
                       data = p2g_data, 
                       family = poisson())
      glm_summ <- summary(model_glm)
      
      glm_base <- glm(gene ~ peak, 
                      data = p2g_data, 
                      family = poisson())
      base_summ <- summary(glm_base)
      
      result <- list(
        gene = o$gene[i],
        peak = o$peak[i],
        p_wald = p_wald, 
        stat_wald = coef_est,
        p_lrt = p_lrt,
        p_glm_basic = base_summ$coefficients["peak", "Pr(>|z|)"],
        stat_glm_basic = base_summ$coefficients["peak", "Estimate"],
        p_peak_glm = glm_summ$coefficients["peak", "Pr(>|z|)"], 
        stat_peak_glm = glm_summ$coefficients["peak", "Estimate"],
        p_condition_glm = glm_summ$coefficients["conditionAct_IL1B_IL23_PGE2", "Pr(>|z|)"], 
        stat_condition_glm = glm_summ$coefficients["conditionAct_IL1B_IL23_PGE2", "Estimate"],
        p_interaction_glm = if("peak:conditionAct_IL1B_IL23_PGE2" %in% rownames(glm_summ$coefficients)) {
          glm_summ$coefficients["peak:conditionAct_IL1B_IL23_PGE2", "Pr(>|z|)"]
        } else {
          NA
        }, 
        stat_interaction_glm = if("peak:conditionAct_IL1B_IL23_PGE2" %in% rownames(glm_summ$coefficients)) {
          glm_summ$coefficients["peak:conditionAct_IL1B_IL23_PGE2", "Estimate"]
        } else {
          NA
        })
      print(i)
      return(result)
    }
  })
  
  results_df <- do.call(rbind, results)
  
  return(results)
}


tcell_dogma <- readRDS(paste0(DATA_PATH, "tcell_annotated_updated.RDS"))
tcell_dogma <- subset(tcell_dogma, subset = condition %in% c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2"))
# tcell_dogma@assays <- tcell_dogma@assays[c("RNA", "peaks")]
DefaultAssay(tcell_dogma) <- "RNA"
gc()

# tcell_dogma <- subset(tcell_dogma, cells = sample(colnames(tcell_dogma), 1000))

# Select Genes
top_gene <- read_xlsx(paste0(WORK_PATH, "Colocalization.xlsx"), 1)
gene_list <- strsplit(top_gene$gene, "\\|")
# gene_vector <- unique(unlist(gene_list))
gene_vector <- c("PTGER4", "CXCR5")
tcell_dogma <- FindVariableFeatures(object = tcell_dogma)
gc()

# no_cores <- 60
# p2g_results_list <- mclapply(unique(tcell_dogma$celltype_updated), function(celltype) {
# 
#   tcell_dogma_subset <- subset(tcell_dogma, subset = celltype_updated == celltype)
#   tcell_dogma_subset <- FindVariableFeatures(tcell_dogma_subset)
#   rna.mat <- tcell_dogma_subset@assays$RNA@data[gene_vector,]
#   peak.mat <- tcell_dogma_subset@assays$peaks@data
#   condition.mat <- data.frame(Barcode = rownames(tcell_dogma_subset@meta.data),
#                               Condition = tcell_dogma_subset$condition)
# 
#   p2g_glmm <- peak2gene_glmm(peak.mat = peak.mat,
#                              rna.mat = rna.mat,
#                              max.dist = 50000, # 50kb
#                              condition.mat = condition.mat)
#   p2g_glmm$celltype <- celltype
#   saveRDS(p2g_glmm, paste0(WORK_PATH, "P2G_", celltype, "_GLMM.rds"))
# 
#   return(p2g_glmm)
# }, mc.cores = no_cores)

# p2g_result <- do.call(rbind, p2g_results_list)

p2g_result <- data.frame()
for(c in c("CD4+ Memory (Activated) - Th17", "CD4+ Memory (Resting) - Th17",
           "CD4+ Memory (Activated) - Th1", "CD4+ Memory (Resting) - Th1")){
  
  tcell_dogma_subset <- subset(tcell_dogma, subset = celltype_updated == c)
  tcell_dogma_subset <- FindVariableFeatures(tcell_dogma_subset)
  rna.mat <- tcell_dogma_subset@assays$RNA@data[gene_vector,]
  peak.mat <- tcell_dogma_subset@assays$peaks@data
  condition.mat <- data.frame(Barcode = rownames(tcell_dogma_subset@meta.data),
                              Condition = tcell_dogma_subset$condition)
  
  p2g_glmm <- peak2gene_glmm(peak.mat = peak.mat,
                             rna.mat = rna.mat,
                             max.dist = 50000, # 50kb
                             condition.mat = condition.mat)
  
  list_to_df <- function(item) {
    as.data.frame(t(unlist(item)))  
  }
  p2g_glmm <- Filter(function(x) !is.null(x), p2g_glmm)
  p2g_glmm <- do.call(rbind, lapply(p2g_glmm, list_to_df))
  p2g_glmm$celltype <- c
  # saveRDS(p2g_glmm, paste0(WORK_PATH, "P2G_", celltype, "_GLMM.rds"))
  p2g_result <- rbind(p2g_result, p2g_glmm)
}

saveRDS(p2g_result, paste0(WORK_PATH, "P2G_GLMM.rds"))


list_to_df <- function(item) {
  as.data.frame(t(unlist(item)))  
}
results_df <- do.call(rbind, lapply(p2g_glmm, list_to_df))
rownames(results_df) <- NULL
