rm(list=ls())
gc()

# Load Packages
library(Seurat)
library(Signac)
library(ggplot2)
library(tidyverse)
library(pROC)
library(caret)

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
        auc_value <- as.numeric(auc_obj)
        
        # Calculate confusion matrix and retrieve PPV and FPR
        confusion <- confusionMatrix(as.factor(class_y), as.factor(class_x))
        ppv <- confusion$byClass['Pos Pred Value']
        specificity <- confusion$byClass['Specificity']
        sensitivity <- confusion$byClass['Sensitivity']
        
      } else {
        auc_value <- NA
        ppv <- NA
        specificity <- NA  # Set FPR to NA if not computable
        sensitivity <- NA
      }
      
      data.frame(
        CellType = unique(.$cell_type),
        SpearmanCorrelation = cor_spearman,
        AUC = auc_value,
        PPV = ppv,
        Sensitivity = sensitivity,
        Specificity = specificity,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        CellType = unique(.$cell_type),
        SpearmanCorrelation = NA,
        AUC = NA,
        PPV = NA,
        FPR = NA,
        stringsAsFactors = FALSE
      )
    }
  }) %>%
  ungroup()
write.csv(results_by_celltype, paste0(WORK_PATH, "Predicted_P2G.csv"))

Results <- Results %>% 
  mutate(class.x = ifelse(EmpFDR.x < 0.05, 1, 0),
         class.y = ifelse(EmpFDR.y < 0.05, 1, 0))

roc_list <- list()
color_palette <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", 
                   "#91D1C2FF", "#DC0000FF", "#7E6148FF",
                   "#3B499299", "#EE000099", "#008B4599", "#63187999", "#00828099", 
                   "#BB002199", "#A2005699", "#1B191999",
                   "#E64B3599", "#FF7F0E99", "#EFC00099")
auc_values <- numeric(length(unique(Results$cell_type)))

for(i in seq_along(unique(Results$cell_type))) {
  cell_type <- unique(Results$cell_type)[i]
  subset_data <- Results[Results$cell_type == cell_type, ]
  roc_obj <- roc(subset_data$class.x, subset_data$class.y)
  auc_value <- auc(roc_obj)
  roc_list[[cell_type]] <- list(roc = roc_obj, auc = auc_value)
  auc_values[i] <- auc_value
}

pdf(paste0(WORK_PATH, "Predicted_P2G.pdf"))
plot(roc_list[[1]]$roc, main="ROC Curve for Peak-to-gene Linkage from Predicted RNA", 
     col=color_palette[1], xlim=c(0, 1), ylim=c(0, 1))

for(j in 2:length(roc_list)) {
  plot(roc_list[[names(roc_list)[j]]]$roc, add=TRUE, col=color_palette[j])
}

auc_min <- min(auc_values, na.rm = TRUE)
auc_max <- max(auc_values, na.rm = TRUE)
text(x=0.2, y=0.1, paste("AUC:", round(auc_min, 3), "to", round(auc_max, 3)), 
     cex=1.2, col="black")
dev.off()


# P2G vs Shuffle
Results <- inner_join(P2G[,c(4:6, 13, 14)], Shuffle[,c(4:6, 13, 14)], 
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
        auc_value <- as.numeric(auc_obj)
        
        # Calculate confusion matrix and retrieve PPV and FPR
        confusion <- confusionMatrix(as.factor(class_y), as.factor(class_x))
        ppv <- confusion$byClass['Pos Pred Value']
        specificity <- confusion$byClass['Specificity']
        sensitivity <- confusion$byClass['Sensitivity']
        
      } else {
        auc_value <- NA
        ppv <- NA
        specificity <- NA  # Set FPR to NA if not computable
        sensitivity <- NA
      }
      
      data.frame(
        CellType = unique(.$cell_type),
        SpearmanCorrelation = cor_spearman,
        AUC = auc_value,
        PPV = ppv,
        Sensitivity = sensitivity,
        Specificity = specificity,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        CellType = unique(.$cell_type),
        SpearmanCorrelation = NA,
        AUC = NA,
        PPV = NA,
        FPR = NA,
        stringsAsFactors = FALSE
      )
    }
  }) %>%
  ungroup()
write.csv(results_by_celltype, paste0(WORK_PATH, "Shuffled_P2G.csv"))

Results <- Results %>% 
  mutate(class.x = ifelse(EmpFDR.x < 0.05, 1, 0),
         class.y = ifelse(EmpFDR.y < 0.05, 1, 0))

roc_list <- list()
color_palette <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", 
                   "#91D1C2FF", "#DC0000FF", "#7E6148FF",
                   "#3B499299", "#EE000099", "#008B4599", "#63187999", "#00828099", 
                   "#BB002199", "#A2005699", "#1B191999",
                   "#E64B3599", "#FF7F0E99", "#EFC00099")
auc_values <- numeric(length(unique(Results$cell_type)))

for(i in seq_along(unique(Results$cell_type))) {
  cell_type <- unique(Results$cell_type)[i]
  subset_data <- Results[Results$cell_type == cell_type, ]
  roc_obj <- roc(subset_data$class.x, subset_data$class.y)
  auc_value <- auc(roc_obj)
  roc_list[[cell_type]] <- list(roc = roc_obj, auc = auc_value)
  auc_values[i] <- auc_value
}

pdf(paste0(WORK_PATH, "Shuffled_P2G.pdf"))
plot(roc_list[[1]]$roc, main="ROC Curve for Peak-to-gene Linkage from Shuffled RNA", 
     col=color_palette[1], xlim=c(0, 1), ylim=c(0, 1))

for(j in 2:length(roc_list)) {
  plot(roc_list[[names(roc_list)[j]]]$roc, add=TRUE, col=color_palette[j])
}

auc_min <- min(auc_values, na.rm = TRUE)
auc_max <- max(auc_values, na.rm = TRUE)
text(x=0.2, y=0.1, paste("AUC:", round(auc_min, 3), "to", round(auc_max, 3)), 
     cex=1.2, col="black")
dev.off()

# # Correlation
# cor(Results$Correlation.x, Results$Correlation.y, method = "spearman")
# Results <- Results %>% 
#   mutate(class.x = ifelse(EmpFDR.x < 0.05, 1, 0),
#          class.y = ifelse(EmpFDR.y < 0.05, 1, 0))
# 
# roc_obj <- roc(Results$class.x, Results$class.y)
# auc_value <- auc(roc_obj)
# 
# plot(roc_obj, main="ROC Curve for Peak-to-gene Linkage from Shuffled RNA", col="#1c61b6")
# auc(roc_obj)
# text(x=0.6, y=0.4, paste("AUC:", round(auc(roc_obj), 3)), cex=1.2)

