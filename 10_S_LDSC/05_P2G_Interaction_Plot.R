rm(list=ls())
gc()

# Load Packages
library(Seurat)
library(Signac)
library(ggplot2)
library(tidyverse)

# Parameters
DATA_PATH <- "/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/04_Global_Analysis/03_P2G/"

# Function
list_to_df <- function(item) {
  as.data.frame(t(unlist(item)))  
}

# Load the results
# Typo: intersection -> interaction
Th1_A <- readRDS(paste0(WORK_PATH, "P2G_CD4+ Memory (Activated) - Th1_GLMM.rds"))
Th1_A[6500] <- NULL
Th1_A <- do.call(rbind, lapply(Th1_A, list_to_df))
Th1_A$celltype <- "Th1_Activated"
Th1_A <- Th1_A %>%
  mutate(
    FDR_wald = p.adjust(p_wald, method = "fdr"),
    FDR_lrt = p.adjust(p_lrt, method = "fdr"),
    FDR_glm_basic = p.adjust(p_glm_basic, method = "fdr"),
    FDR_peak_glm = p.adjust(p_peak_glm, method = "fdr"),
    FDR_condition_glm = p.adjust(p_condition_glm, method = "fdr"),
    FDR_intersection_glm = p.adjust(p_intersection_glm, method = "fdr")
  )

Th1_R <- readRDS(paste0(WORK_PATH, "P2G_CD4+ Memory (Resting) - Th1_GLMM.rds"))
Th1_R[6500] <- NULL
Th1_R <- do.call(rbind, lapply(Th1_R, list_to_df))
Th1_R$celltype <- "Th1_Resting"
Th1_R<- Th1_R %>%
  mutate(
    FDR_wald = p.adjust(p_wald, method = "fdr"),
    FDR_lrt = p.adjust(p_lrt, method = "fdr"),
    FDR_glm_basic = p.adjust(p_glm_basic, method = "fdr"),
    FDR_peak_glm = p.adjust(p_peak_glm, method = "fdr"),
    FDR_condition_glm = p.adjust(p_condition_glm, method = "fdr"),
    FDR_intersection_glm = p.adjust(p_intersection_glm, method = "fdr")
  )

Th17_R <- readRDS(paste0(WORK_PATH, "P2G_CD4+ Memory (Resting) - Th17_GLMM.rds"))
Th17_R[6500] <- NULL
Th17_R <- do.call(rbind, lapply(Th17_R, list_to_df))
Th17_R$celltype <- "Th17_Resting"
Th17_R <- Th17_R %>%
  mutate(
    FDR_wald = p.adjust(p_wald, method = "fdr"),
    FDR_lrt = p.adjust(p_lrt, method = "fdr"),
    FDR_glm_basic = p.adjust(p_glm_basic, method = "fdr"),
    FDR_peak_glm = p.adjust(p_peak_glm, method = "fdr"),
    FDR_condition_glm = p.adjust(p_condition_glm, method = "fdr"),
    FDR_intersection_glm = p.adjust(p_intersection_glm, method = "fdr")
  )

Th17_A <- readRDS(paste0(WORK_PATH, "P2G_CD4+ Memory (Activated) - Th17_GLMM.rds"))
Th17_A[6500] <- NULL
Th17_A <- do.call(rbind, lapply(Th17_A, list_to_df))
Th17_A$celltype <- "Th17_Activated"
Th17_A <- Th17_A %>%
  mutate(
    FDR_wald = p.adjust(p_wald, method = "fdr"),
    FDR_lrt = p.adjust(p_lrt, method = "fdr"),
    FDR_glm_basic = p.adjust(p_glm_basic, method = "fdr"),
    FDR_peak_glm = p.adjust(p_peak_glm, method = "fdr"),
    FDR_condition_glm = p.adjust(p_condition_glm, method = "fdr"),
    FDR_intersection_glm = p.adjust(p_intersection_glm, method = "fdr")
  )

# Combine
Results <- rbind(Th1_A, Th1_R, Th17_A, Th17_R)

# Summary
## 25995 peak-gene pairs
sum(Results$FDR_wald[!is.na(Results$FDR_wald)] < 0.05) # 24991 GLMM
sum(Results$FDR_lrt[!is.na(Results$FDR_lrt)] < 0.05) # 24412 GLMM
sum(Results$FDR_glm_basic < 0.05) # 25448 GLM Basic
sum(Results$FDR_peak_glm < 0.05) # 24936 GLM Peak
sum(Results$FDR_condition_glm < 0.05) # 14568 GLM Condition
sum(Results$FDR_intersection_glm < 0.05) # 12268 GLM Interaction

Results_sig <- Results %>% filter(FDR_glm_basic < 0.05)
sum(Results_sig$FDR_peak_glm < 0.05) # GLM Peak
sum(Results_sig$FDR_condition_glm < 0.05) # GLM Condition
sum(Results_sig$FDR_intersection_glm < 0.05) # GLM Interaction

# Load data
tcell_dogma <- readRDS(paste0(DATA_PATH, "tcell_annotated_updated.RDS"))
tcell_dogma <- subset(tcell_dogma, subset = condition %in% c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2"))
tcell_dogma <- subset(tcell_dogma, subset = celltype_updated %in% c("CD4+ Memory (Activated) - Th17", 
                                                                    "CD4+ Memory (Resting) - Th17", 
                                                                    "CD4+ Memory (Activated) - Th1", 
                                                                    "CD4+ Memory (Resting) - Th1"))
condition.mat <- data.frame(Barcode = rownames(tcell_dogma@meta.data),
                            CellType = tcell_dogma$celltype_updated,
                            Condition = tcell_dogma$condition)

peak <- "chr5-40652340-40653946"
gene <- "PTGER4"
a <- Results[which(Results$gene == gene & Results$peak == peak),]
rna.mat <- tcell_dogma@assays$RNA@data[gene,]
peak.mat <- tcell_dogma@assays$peaks@data[peak,]

# GLM
barcode <- condition.mat$Barcode[condition.mat$CellType == "CD4+ Memory (Activated) - Th1"]
p2g_data <- data.frame(gene = rna.mat[barcode],
                       peak = peak.mat[barcode],
                       condition = condition.mat$Condition[condition.mat$Barcode %in% barcode])
model_glm <- glm(gene ~ peak * condition, 
                 data = p2g_data, 
                 family = poisson())
summary(model_glm)

# Start from here
p2g_result <- p2g_result %>%
  mutate(
    FDR_wald = p.adjust(p_wald, method = "fdr"),
    FDR_lrt = p.adjust(p_lrt, method = "fdr"),
    FDR_glm_basic = p.adjust(p_glm_basic, method = "fdr"),
    FDR_peak_glm = p.adjust(p_peak_glm, method = "fdr"),
    FDR_condition_glm = p.adjust(p_condition_glm, method = "fdr"),
    FDR_interaction_glm = p.adjust(p_interaction_glm, method = "fdr")
  )
## 156 peak-gene pairs
sum(p2g_result$FDR_wald[!is.na(p2g_result$FDR_wald)] < 0.05) # 97 GLMM
sum(p2g_result$FDR_lrt[!is.na(p2g_result$FDR_lrt)] < 0.05) # 79 GLMM
sum(p2g_result$FDR_glm_basic < 0.05) # 116 GLM Basic
sum(p2g_result$FDR_peak_glm < 0.05) # 72 GLM Peak
sum(p2g_result$FDR_condition_glm < 0.05) # 71 GLM Condition
sum(p2g_result$FDR_interaction_glm[!is.na(p2g_result$FDR_interaction_glm)] < 0.05) # 15 GLM Interaction

a <- p2g_result[p2g_result$FDR_glm_basic < 0.05,]
# Scatter Plot
peak <- c("chr5-40652340-40653946", "chr11-118928141-118930490")
gene <- c("PTGER4", "CXCR5")
rna.mat <- tcell_dogma@assays$RNA@data[gene,]
peak.mat <- tcell_dogma@assays$peaks@data[peak,]
p2g_data <- rbind(t(rna.mat), t(peak.mat)) %>% data.frame()


