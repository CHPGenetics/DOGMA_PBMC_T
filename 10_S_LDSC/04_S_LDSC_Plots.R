rm(list=ls())
gc()

# Load Packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(data.table)

# Parameters
DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Output/"
WORK_PATH <- paste0(DATA_PATH, "Plots/")

# IBD DE
IBD_Seurat_DE <- fread(paste0(DATA_PATH, "Seurat_DE/Result/IBD_S_LDSC_Result.cell_type_results.txt"))
IBD_Seurat_DE$Method <- "Seurat_DE"
IBD_Pseudo_Condition_DE <- fread(paste0(DATA_PATH, "Pseudo_Condition_DE/Result/IBD_S_LDSC_Result.cell_type_results.txt"))
IBD_Pseudo_Condition_DE$Method <- "Pseudo_Condition_DE"


# IBD DA
IBD_Pseudo_Condition_DA <- fread(paste0(DATA_PATH, "Pseudo_Condition_DA/Result/IBD_S_LDSC_Result.cell_type_results.txt"))


# UC DE
UC_Seurat_DE <- fread(paste0(DATA_PATH, "Seurat_DE/Result/UC_S_LDSC_Result.cell_type_results.txt"))
UC_Pseudo_Condition_DE <- fread(paste0(DATA_PATH, "Pseudo_Condition_DE/Result/UC_S_LDSC_Result.cell_type_results.txt"))

# UC DA


# CD DE
CD_Seurat_DE <- fread(paste0(DATA_PATH, "Seurat_DE/Result/CD_S_LDSC_Result.cell_type_results.txt"))
CD_Seurat_DE$Method <- "Seurat_DE"
CD_Pseudo_Condition_DE <- fread(paste0(DATA_PATH, "Pseudo_Condition_DE/Result/CD_S_LDSC_Result.cell_type_results.txt"))
CD_Pseudo_Condition_DE$Method <- "Pseudo_Condition_DE"
combin_output <- rbind(CD_Seurat_DE, CD_Pseudo_Condition_DE)
combin_output$Name <- gsub("_", " ", combin_output$Name)
combin_output$Name <- factor(combin_output$Name, levels = c("CD4+ Naive (Resting)", "CD4+ Naive (Activated)", "CD4+ Regulatory (Resting)",
                                                            "CD4+ Regulatory (Activated)", "CD4+ Memory (Resting) - Th1", "CD4+ Memory (Activated) - Th1",
                                                            "CD4+ Memory (Resting) - Th17", "CD4+ Memory (Activated) - Th17", "CD4+ Memory (Resting) - Tfh",
                                                            "CD4+ Memory (Activated) - Tfh", "CD4+ Memory (Resting) - Other", "CD4+ Memory (Activated) - Other",
                                                            "CD8+ Naive (Resting)", "CD8+ Naive (Activated)", "CD8+ Regulatory",
                                                            "CD8+ Memory (Resting)", "CD8+ Memory (Activated)", "MAITs (Resting)", 
                                                            "MAITs (Activated)", "Gamma Delta"))

ggplot(data = combin_output, aes(x = Name, y = -log10(Coefficient_P_value))) + 
  geom_bar(aes(fill = Method), stat = 'identity', position = position_dodge(0.8), 
           width = 0.8, show.legend = T) +
  scale_fill_manual(values = c("#DC0000B2","#4DBBD5B2")) + 
  geom_hline(yintercept = -log10(0.05/nrow(combin_output)), color = "gold", lty = 5, lwd = 1) + 
  labs(x = "cell types", y = "-log10(P value)") +
  scale_y_continuous(limits = c(0, 6), breaks = c(seq(0, 6, 0.5))) +
  coord_flip() +
  ggtitle("CD") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(size = 15, color = "grey20"),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 18, color = "grey20"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

# CD DA


