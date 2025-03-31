rm(list=ls())
gc()

# Load Packages
library(Seurat)
library(Signac)
library(ggplot2)
library(tidyverse)
library(pROC)
library(ggsci)

# Parameters
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/04_Global_Analysis/04_P2G_Shuffle/"

# Load data
P2G <- readRDS(paste0(WORK_PATH, "P2G/Emprical/All_P2G.rds"))
Predict <- readRDS(paste0(WORK_PATH, "P2G_Predicted/Emprical/All_P2G.rds"))
Shuffle <- readRDS(paste0(WORK_PATH, "P2G_Shuffle/Emprical/All_P2G.rds"))

# P2G vs Predicted
Results <- inner_join(P2G[,c(4:6, 13, 14)], Predict[,c(4:6, 13, 14)], 
                      by = c("gene", "peak", "cell_type"))
Results <- Results %>%
  filter(Correlation.x != "NaN" & Correlation.y != "NaN")

# Correlation
results_by_celltype <- Results %>%
  group_by(cell_type) %>%
  do({
    if(nrow(.) >= 2) {
      cor_spearman <- cor(.$Correlation.x, .$Correlation.y, method = "spearman")
      
      class_x <- ifelse(.$EmpFDR.x < 0.05, 1, 0)
      class_y <- ifelse(.$EmpFDR.y < 0.05, 1, 0)
      
      if(length(class_x) == length(class_y)) {
        roc_obj <- roc(class_x, class_y)
        auc_obj <- auc(roc_obj)
        auc_value <- as.numeric(auc_obj)  # Convert AUC to numeric if not already
      } else {
        auc_value <- NA
      }
      
      data.frame(
        CellType = unique(.$cell_type),
        SpearmanCorrelation = cor_spearman,
        AUC = auc_value,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        CellType = unique(.$cell_type),
        SpearmanCorrelation = NA,
        AUC = NA,
        stringsAsFactors = FALSE
      )
    }
  }) %>%
  ungroup()

cor(Results$Correlation.x, Results$Correlation.y, method = "spearman")
Results <- Results %>% 
  mutate(class.x = ifelse(EmpFDR.x < 0.05, 1, 0),
         class.y = ifelse(EmpFDR.y < 0.05, 1, 0))

roc_list <- list()
colors <- rainbow(length(unique(Results$cell_type)))
for(i in unique(Results$cell_type)) {
  subset_data <- Results[Results$cell_type == i, ]
  roc_obj <- roc(subset_data$class.x, subset_data$class.y)
  auc_value <- auc(roc_obj)
  roc_list[[i]] <- list(roc = roc_obj, auc = auc_value)
}

plot(roc_list[[1]]$roc, main="ROC Curve for Peak-to-gene Linkage from Predicted RNA", col=colors[1], xlim=c(0, 1), ylim=c(0, 1))
text(x=0.2, y=0.4, paste("AUC:", round(roc_list[[1]]$auc, 3)), cex=0.8, col=colors[1])

for(j in 2:length(roc_list)) {
  plot(roc_list[[names(roc_list)[j]]]$roc, add=TRUE, col=colors[j])
  text(x=0.2, y=0.4 - 0.05 * j, paste("AUC:", round(roc_list[[names(roc_list)[j]]]$auc, 3)), cex=0.8, col=colors[j])
}




roc_obj <- roc(Results$class.x, Results$class.y)
auc_value <- auc(roc_obj)

plot(roc_obj, main="ROC Curve for Peak-to-gene Linkage from Predicted RNA", col="#1c61b6")
auc(roc_obj)
text(x=0.6, y=0.4, paste("AUC:", round(auc(roc_obj), 3)), cex=1.2)


# P2G vs Shuffle
Results <- inner_join(P2G[,c(4:6, 13, 14)], Shuffle[,c(4:6, 13, 14)], 
                      by = c("gene", "peak", "cell_type"))
Results <- Results %>%
  filter(Correlation.x != "NaN" & Correlation.y != "NaN")

# Correlation
cor(Results$Correlation.x, Results$Correlation.y, method = "spearman")
Results <- Results %>% 
  mutate(class.x = ifelse(EmpFDR.x < 0.05, 1, 0),
         class.y = ifelse(EmpFDR.y < 0.05, 1, 0))

roc_obj <- roc(Results$class.x, Results$class.y)
auc_value <- auc(roc_obj)

plot(roc_obj, main="ROC Curve for Peak-to-gene Linkage from Shuffled RNA", col="#1c61b6")
auc(roc_obj)
text(x=0.6, y=0.4, paste("AUC:", round(auc(roc_obj), 3)), cex=1.2)

confusion_matrix <- table(Actual = Results$class.x, Predicted = Results$class.y)
true_positives <- confusion_matrix[2, 2]
false_negatives <- confusion_matrix[2, 1]
true_positive_rate <- true_positives / (true_positives + false_negatives)
false_negative_rate <- false_negatives / (true_positives + false_negatives)

