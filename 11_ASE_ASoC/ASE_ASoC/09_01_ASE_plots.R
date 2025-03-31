rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(UpSetR)
library(ggvenn)

# Parameters
DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASE_Count/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/"

# Load data
ASE_Results <- fread(paste0(DATA_PATH, "ASE_counts.txt"))
ASE_Results <- ASE_Results %>%
  filter(p_value_adj < 0.05) %>%
  mutate(sample_id = paste0(condition, "_", sample))

# Summary
PGE2_ASE <- ASE_Results %>% filter(condition == "Act_IL1B_IL23_PGE2")
TGFB_ASE <- ASE_Results %>% filter(condition == "Act_IL1B_IL23_PGE2_TGFB")
TWO_ASE <- ASE_Results %>% filter(condition == "Act_IL1B_IL23_TGFB")
Control_ASE <- ASE_Results %>% filter(condition == "Act_IL1B_IL23")

# Proportion
PGE2_unique_snp <- setdiff(PGE2_ASE$variantID, Control_ASE$variantID) %>%
  setdiff(TGFB_ASE$variantID) %>% setdiff(TWO_ASE$variantID)
TGFB_unique_snp <- setdiff(TGFB_ASE$variantID, Control_ASE$variantID) %>%
  setdiff(PGE2_ASE$variantID) %>% setdiff(TWO_ASE$variantID)
Control_unique_snp <- setdiff(Control_ASE$variantID, PGE2_ASE$variantID) %>%
  setdiff(TGFB_ASE$variantID) %>% setdiff(TWO_ASE$variantID)
TWO_unique_snp <- setdiff(TWO_ASE$variantID, PGE2_ASE$variantID) %>%
  setdiff(TGFB_ASE$variantID) %>% setdiff(Control_ASE$variantID)
Shared_snp <- intersect(unique(PGE2_ASE$variantID), unique(Control_ASE$variantID)) %>%
  intersect(unique(TGFB_ASE$variantID)) %>% intersect(unique(TWO_ASE$variantID))


## Top ASE cites
top_ASE <- table(ASE_Results$variantID) %>% 
  data.frame() %>%
  dplyr::arrange(desc(Freq))

## Top ASE cites in Control
top_ASE_Control <- ASE_Results %>%
  filter(variantID %in% Control_unique_snp)
top_ASE_Control <- table(top_ASE_Control$variantID) %>% 
  data.frame() %>%
  dplyr::arrange(desc(Freq))

## Top ASE cites in PGE2
top_ASE_PGE2 <- ASE_Results %>%
  filter(variantID %in% PGE2_unique_snp)
top_ASE_PGE2 <- table(top_ASE_PGE2$variantID) %>% 
  data.frame() %>%
  dplyr::arrange(desc(Freq))

## Top Shared ASE cites
top_ASE_shared <- ASE_Results %>%
  filter(variantID %in% Shared_snp)
top_ASE_shared <- table(top_ASE_shared$variantID) %>% 
  data.frame() %>%
  dplyr::arrange(desc(Freq))

# Forest plot
ASE_Results <- ASE_Results %>%
  mutate(ref_prop = refCount/totalCount)

## Shared
ASE_Results_sub <- ASE_Results %>%
  filter(variantID %in% Shared_snp) %>%
  filter(variantID %in% top_ASE$Var1[1:21]) %>%
  mutate(variantID = factor(variantID, levels = rev(top_ASE$Var1))) %>%
  arrange(variantID)

pdf(paste0(WORK_PATH, "ASE_Ref_proportion_shared.pdf"), width = 6, height = 5)
ggplot(ASE_Results_sub, aes(y = variantID, x = ref_prop, color = condition)) +
  geom_point() +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  scale_color_manual(values = c("#4DBBD5B2", "#DC0000B2", "#a19ba9", "#7a929e"),
                     labels = c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2",
                                "Act_IL1B_IL23_PGE2_TGFB", "Act_IL1B_IL23_TGFB")) +
  xlim(0, 1) +
  theme_bw() +
  labs(x = "Reference Allele Mapping Proportion", y = "SNPs")
dev.off()

## Unique (PGE2)
ASE_Results_sub <- ASE_Results %>%
  filter(variantID %in% PGE2_unique_snp) %>%
  filter(variantID %in% top_ASE_PGE2$Var1[1:10]) %>%
  mutate(variantID = factor(variantID, levels = rev(top_ASE_PGE2$Var1))) %>%
  arrange(variantID)

pdf(paste0(WORK_PATH, "ASE_Ref_proportion_PGE2.pdf"), width = 6, height = 5)
ggplot(ASE_Results_sub, aes(y = variantID, x = ref_prop, color = condition)) +
  geom_point() +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  scale_color_manual(values = c("#DC0000B2"),
                     labels = c("Act_IL1B_IL23_PGE2")) +
  xlim(0, 1) +
  theme_bw() +
  labs(x = "Reference Allele Mapping Proportion", y = "SNPs")
dev.off()

## Unique (Control)
ASE_Results_sub <- ASE_Results %>%
  filter(variantID %in% Control_unique_snp) %>%
  filter(variantID %in% top_ASE_Control$Var1[1:10]) %>%
  mutate(variantID = factor(variantID, levels = rev(top_ASE_Control$Var1))) %>%
  arrange(variantID)

pdf(paste0(WORK_PATH, "ASE_Ref_proportion_Control.pdf"), width = 6, height = 5)
ggplot(ASE_Results_sub, aes(y = variantID, x = ref_prop, color = condition)) +
  geom_point() +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  scale_color_manual(values = c("#4DBBD5B2"),
                     labels = c("Act_IL1B_IL23")) +
  xlim(0, 1) +
  theme_bw() +
  labs(x = "Reference Allele Mapping Proportion", y = "SNPs")
dev.off()

# Upset plot
input <- c(Control = length(unique(ASE_Results$variantID[which(ASE_Results$condition == "Act_IL1B_IL23")])),
           PGE2 = length(unique(ASE_Results$variantID[which(ASE_Results$condition == "Act_IL1B_IL23_")])),
           TGFB = length(unique(ASE_Results$variantID[which(ASE_Results$condition == "Act_IL1B_IL23_TGFB")])),
           PGE2_TGFB = length(unique(ASE_Results$variantID[which(ASE_Results$condition == "Act_IL1B_IL23_PGE2_TGFB")])),
           "Control&PGE2" = length(intersect(unique(PGE2_ASE$variantID), unique(Control_ASE$variantID))),
           "Control&TGFB" = length(intersect(unique(TGFB_ASE$variantID), unique(Control_ASE$variantID))),
           "Control&PGE2_TGFB" = length(intersect(unique(TWO_ASE$variantID), unique(Control_ASE$variantID))),
           "PGE2&TGFB" = length(intersect(unique(PGE2_ASE$variantID), unique(TGFB_ASE$variantID))),
           "PGE2&PGE2_TGFB" = length(intersect(unique(PGE2_ASE$variantID), unique(TWO_ASE$variantID))),
           "TGFB&PGE2_TGFB" = length(intersect(unique(TGFB_ASE$variantID), unique(TWO_ASE$variantID))),
           "Control&PGE2&TGFB&PGE2_TGFB" = length(Shared_snp))

pdf(paste0(WORK_PATH, "ASE_Upset_plot.pdf"), width = 12, height = 6)
upset(fromExpression(input), 
      nintersects = 40, 
      nsets = 10, 
      order.by = "freq", 
      decreasing = TRUE, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 2, 
      point.size = 3.8, 
      line.size = 1.5)
dev.off()

# Venn plot
x <- list(Act_IL1B_IL23_PGE2 = unique(ASE_Results$variantID[which(ASE_Results$condition == "Act_IL1B_IL23_PGE2")]),
          Act_IL1B_IL23_PGE2_TGFB = unique(ASE_Results$variantID[which(ASE_Results$condition == "Act_IL1B_IL23_PGE2_TGFB")]),
          Act_IL1B_IL23_TGFB = unique(ASE_Results$variantID[which(ASE_Results$condition == "Act_IL1B_IL23_TGFB")]),
          Act_IL1B_IL23 = unique(ASE_Results$variantID[which(ASE_Results$condition == "Act_IL1B_IL23")]))

pdf(paste0(WORK_PATH, "ASE_Venn_plot.pdf"), width = 11, height = 6.8)
ggvenn(x, 
       fill_color = c("#CD534CFF", "#0073C2FF", "#EFC000FF", "#868686FF"),
       stroke_size = 0.5, set_name_size = 3) # Reduce set label size
dev.off()

###############################################################################################
# Cell-type specific
## Overall summary Bar chart
ASE_Results <- fread(paste0("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASE_Count/ASE_counts.txt")) %>%
  filter(p_value_adj < 0.05)
ASoC_Results <- fread(paste0("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASoC_Count/ASoC_counts.txt")) %>%
  filter(p_value_adj < 0.05)

library(openxlsx)
write.xlsx(ASE_Results, file = "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASE_Count/ASE.xlsx")
write.xlsx(ASoC_Results, file = "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASE_Count/ASoC.xlsx")

ASE_ASoC <- data.frame()
for(i in unique(ASoC_Results$celltype)){
  
  for(j in unique(ASoC_Results$condition)){
    
    ASE_subset <- ASE_Results %>% filter(celltype == i) %>%
      filter(condition == j)
    ASE_ASoC <- rbind(ASE_ASoC, c(i, length(unique(ASE_subset$variantID)), "ASE", j))
  }
}
for(i in unique(ASoC_Results$celltype)){
  
  for(j in unique(ASoC_Results$condition)){
    
    ASoC_subset <- ASoC_Results %>% filter(celltype == i) %>%
      filter(condition == j)
    ASE_ASoC <- rbind(ASE_ASoC, c(i, length(unique(ASoC_subset$variantID)), "ASoC", j))
  }
}
colnames(ASE_ASoC) <- c("cell_type", "value", "measurement", "condition")
ASE_ASoC$cell_type <- gsub("_", " ", ASE_ASoC$cell_type)
ASE_ASoC$value <- as.numeric(ASE_ASoC$value)
ASE_ASoC$cell_type <- factor(ASE_ASoC$cell_type, levels = c("CD4+ Naive Resting", "CD4+ Naive Activated", 
                                                            "CD4+ Regulatory Resting", "CD4+ Regulatory Activated",
                                                            "CD4+ Memory Resting Th1", "CD4+ Memory Activated Th1",
                                                            "CD4+ Memory Resting Th17", "CD4+ Memory Activated Th17",
                                                            "CD4+ Memory Resting Tfh", "CD4+ Memory Activated Tfh",
                                                            "CD4+ Memory Resting Other", "CD4+ Memory Activated Other",
                                                            "CD8+ Naive Resting", "CD8+ Naive Activated", 
                                                            "CD8+ Regulatory", "CD8+ Memory Resting", 
                                                            "CD8+ Memory Activated", "MAITs Resting", 
                                                            "MAITs Activated", "Gamma Delta"))

data1 <- ASE_ASoC %>% filter(condition == 'Act_IL1B_IL23') %>% filter(measurement == "ASE") %>%
  filter(cell_type %in% c("CD4+ Naive Resting", "CD4+ Naive Activated", 
                          "CD4+ Regulatory Resting", "CD4+ Regulatory Activated",
                          "CD4+ Memory Resting Th1", "CD4+ Memory Activated Th1",
                          "CD4+ Memory Resting Th17", "CD4+ Memory Activated Th17",
                          "CD4+ Memory Resting Tfh", "CD4+ Memory Activated Tfh",
                          "CD4+ Memory Resting Other", "CD4+ Memory Activated Other",
                          "CD8+ Naive Resting", "CD8+ Naive Activated", 
                          "CD8+ Memory Resting", "CD8+ Memory Activated", "MAITs Resting", "MAITs Activated")) %>%
  arrange(cell_type)

wilcox.test(data1$value[seq(from = 2, to = 18, by = 2)], data1$value[seq(from = 1, to = 17, by = 2)],
            alternative = "greater",)

data2 <- ASE_ASoC %>% filter(condition == 'Act_IL1B_IL23') %>% filter(measurement == "ASoC") %>%
  filter(cell_type %in% c("CD4+ Naive Resting", "CD4+ Naive Activated", 
                          "CD4+ Regulatory Resting", "CD4+ Regulatory Activated",
                          "CD4+ Memory Resting Th1", "CD4+ Memory Activated Th1",
                          "CD4+ Memory Resting Th17", "CD4+ Memory Activated Th17",
                          "CD4+ Memory Resting Tfh", "CD4+ Memory Activated Tfh",
                          "CD4+ Memory Resting Other", "CD4+ Memory Activated Other",
                          "CD8+ Naive Resting", "CD8+ Naive Activated", 
                          "CD8+ Memory Resting", "CD8+ Memory Activated", "MAITs Resting", "MAITs Activated")) %>%
  arrange(cell_type)

wilcox.test(data2$value[seq(from = 2, to = 18, by = 2)], data2$value[seq(from = 1, to = 17, by = 2)], 
            alternative = "greater", paired = TRUE)

data1 <- ASE_ASoC %>% filter(condition == 'Act_IL1B_IL23_PGE2') %>% filter(measurement == "ASE") %>%
  filter(cell_type %in% c("CD4+ Naive Resting", "CD4+ Naive Activated", 
                          "CD4+ Regulatory Resting", "CD4+ Regulatory Activated",
                          "CD4+ Memory Resting Th1", "CD4+ Memory Activated Th1",
                          "CD4+ Memory Resting Th17", "CD4+ Memory Activated Th17",
                          "CD4+ Memory Resting Tfh", "CD4+ Memory Activated Tfh",
                          "CD4+ Memory Resting Other", "CD4+ Memory Activated Other",
                          "CD8+ Naive Resting", "CD8+ Naive Activated", 
                          "CD8+ Memory Resting", "CD8+ Memory Activated", "MAITs Resting", "MAITs Activated")) %>%
  arrange(cell_type)

wilcox.test(data1$value[seq(from = 2, to = 18, by = 2)], data1$value[seq(from = 1, to = 17, by = 2)], 
            alternative = "less", paired = TRUE)

data2 <- ASE_ASoC %>% filter(condition == 'Act_IL1B_IL23_PGE2') %>% filter(measurement == "ASoC") %>%
  filter(cell_type %in% c("CD4+ Naive Resting", "CD4+ Naive Activated", 
                          "CD4+ Regulatory Resting", "CD4+ Regulatory Activated",
                          "CD4+ Memory Resting Th1", "CD4+ Memory Activated Th1",
                          "CD4+ Memory Resting Th17", "CD4+ Memory Activated Th17",
                          "CD4+ Memory Resting Tfh", "CD4+ Memory Activated Tfh",
                          "CD4+ Memory Resting Other", "CD4+ Memory Activated Other",
                          "CD8+ Naive Resting", "CD8+ Naive Activated", 
                          "CD8+ Memory Resting", "CD8+ Memory Activated", "MAITs Resting", "MAITs Activated")) %>%
  arrange(cell_type)

wilcox.test(data2$value[seq(from = 2, to = 18, by = 2)], data2$value[seq(from = 1, to = 17, by = 2)], 
            alternative = "less", paired = TRUE)

data1 <- ASE_ASoC %>% filter(condition %in% c('Act_IL1B_IL23_PGE2', 'Act_IL1B_IL23')) %>% 
  filter(measurement == "ASE") %>%
  filter(cell_type %in% c("CD4+ Naive Resting", "CD4+ Naive Activated", 
                          "CD4+ Regulatory Resting", "CD4+ Regulatory Activated",
                          "CD4+ Memory Resting Th1", "CD4+ Memory Activated Th1",
                          "CD4+ Memory Resting Th17", "CD4+ Memory Activated Th17",
                          "CD4+ Memory Resting Tfh", "CD4+ Memory Activated Tfh",
                          "CD4+ Memory Resting Other", "CD4+ Memory Activated Other",
                          "CD8+ Naive Resting", "CD8+ Naive Activated", 
                          "CD8+ Memory Resting", "CD8+ Memory Activated", "MAITs Resting", "MAITs Activated")) %>%
  arrange(cell_type)

# resting
wilcox.test(data1$value[seq(from = 1, to = 33, by = 4)], data1$value[seq(from = 2, to = 34, by = 4)], 
            alternative = "less", paired = TRUE)
# activated
wilcox.test(data1$value[seq(from = 3, to = 35, by = 4)], data1$value[seq(from = 4, to = 36, by = 4)], 
            alternative = "greater", paired = TRUE)

data1 <- ASE_ASoC %>% filter(condition %in% c('Act_IL1B_IL23_PGE2', 'Act_IL1B_IL23')) %>% 
  filter(measurement == "ASoC") %>%
  filter(cell_type %in% c("CD4+ Naive Resting", "CD4+ Naive Activated", 
                          "CD4+ Regulatory Resting", "CD4+ Regulatory Activated",
                          "CD4+ Memory Resting Th1", "CD4+ Memory Activated Th1",
                          "CD4+ Memory Resting Th17", "CD4+ Memory Activated Th17",
                          "CD4+ Memory Resting Tfh", "CD4+ Memory Activated Tfh",
                          "CD4+ Memory Resting Other", "CD4+ Memory Activated Other",
                          "CD8+ Naive Resting", "CD8+ Naive Activated", 
                          "CD8+ Memory Resting", "CD8+ Memory Activated", "MAITs Resting", "MAITs Activated")) %>%
  arrange(cell_type)

# resting
wilcox.test(data1$value[seq(from = 1, to = 33, by = 4)], data1$value[seq(from = 2, to = 34, by = 4)], 
            alternative = "less", paired = TRUE)
# activated
wilcox.test(data1$value[seq(from = 3, to = 35, by = 4)], data1$value[seq(from = 4, to = 36, by = 4)], 
            alternative = "greater", paired = TRUE)


###################################################################################################
ASE_ASoC <- ASE_ASoC %>%
  mutate(text = ifelse(measurement == "ASoC", paste0("          ", value),
                       ifelse(measurement == "ASE", paste0(value, "          "), NA))) %>%
  mutate(value = ifelse(measurement == "ASoC", value,
                        ifelse(measurement == "ASE", -value, NA)))
pdf(paste0(WORK_PATH, "ASE_ASoC_count.pdf"), width = 9, height = 7)
ggplot(ASE_ASoC, aes(x = cell_type, y = value, fill = measurement)) +
  geom_bar(stat = "identity", position = "stack") +  # Stack the bars
  geom_text(aes(label = text, vjust = 0.7), color = "black", size = 3) + # Add text labels
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("ASE" = "#4DBBD5B2", "ASoC" = "#DC0000B2")) +
  labs(y = "Number of allelic imbalance sites", x = "") +
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.border = element_blank(),   
        strip.background = element_blank(),  
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_blank(),          
        axis.ticks.x = element_blank()) +
  scale_y_continuous(labels = function(x) format(abs(x), scientific = FALSE), 
                     limits = c(-6000, 9000)) +
  facet_wrap(~ condition)
dev.off()


###################################################################################################
# Forest plot
variants_select <- c("5_159209504", "6_31356323", "1_16514093", "5_132482293", "6_31354138",
                     "1_16514084", "3_93470628", "5_142109242", "5_132482268", "6_31354249")
ASE_Results <- fread(paste0("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASE_Count/ASE_counts.txt")) %>%
  mutate(measurement = "ASE") %>% filter(variantID %in% variants_select)
ASoC_Results <- fread(paste0("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASoC_Count/ASoC_counts.txt")) %>%
  mutate(measurement = "ASoC") %>% filter(variantID %in% variants_select)
ASE_ASoC <- rbind(ASE_Results, ASoC_Results)

ASE_ASoC <- ASE_ASoC %>%
  mutate(ref_prop = refCount/totalCount)
data_summary <- ASE_ASoC %>%
  group_by(variantID, celltype, measurement, condition) %>%
  summarise(
    mean_prop = mean(ref_prop, na.rm = TRUE),
    n = n(),
    low = min(ref_prop, na.rm = TRUE),  
    high = max(ref_prop, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(n > 2)
data_summary_sub <- data_summary %>%
  filter(variantID %in% c("5_159209504", "1_16514084"))

data_summary_sub$celltype <- gsub("_", " ", data_summary_sub$celltype)
data_summary_sub$celltype <- factor(data_summary_sub$celltype, levels = c("CD4+ Naive Resting", "CD4+ Naive Activated", 
                                                            "CD4+ Regulatory Resting", "CD4+ Regulatory Activated",
                                                            "CD4+ Memory Resting Th1", "CD4+ Memory Activated Th1",
                                                            "CD4+ Memory Resting Th17", "CD4+ Memory Activated Th17",
                                                            "CD4+ Memory Resting Tfh", "CD4+ Memory Activated Tfh",
                                                            "CD4+ Memory Resting Other", "CD4+ Memory Activated Other",
                                                            "CD8+ Naive Resting", "CD8+ Naive Activated", 
                                                            "CD8+ Regulatory", "CD8+ Memory Resting", 
                                                            "CD8+ Memory Activated", "MAITs Resting", 
                                                            "MAITs Activated", "Gamma Delta"))
ggplot(data_summary_sub, aes(x = mean_prop, 
                             y = celltype, color = measurement)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(xmin = low, xmax = high), 
                width = 0.3, 
                position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  xlim(0, 1) +
  scale_color_manual(values = c("ASE" = "#4DBBD5B2", "ASoC" = "#DC0000B2")) +
  facet_grid(variantID ~ condition, scales = "free_x") +
  labs(x = NULL, y = NULL, title = "Forest Plot of Conditions by Cell Type") +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 9),
        legend.position = "bottom")
#####
# Check IL23R (1_67138637-67265903), 
# IL1R1 (1_102070390-102179874) and PTGER4 (5_40679915-40746800)
ASE_Results <- fread(paste0("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASE_Count/ASE_counts.txt")) %>%
  mutate(measurement = "ASE")
ASoC_Results <- fread(paste0("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASoC_Count/ASoC_counts.txt")) %>%
  mutate(measurement = "ASoC")
ASE_ASoC <- rbind(ASE_Results, ASoC_Results) %>%
  filter((contig == "chr5" & position >= 40679915 & position <= 40746800) |
           (contig == "chr1" & position >= 102070390 & position <= 102179874) | 
           (contig == "chr1" & position >= 67138637 & position <= 67265903))

variants_select <- c('1_67172321', '5_40685693', '5_40687361', '5_40687957', '5_40692838')
ASE_ASoC <- ASE_ASoC %>% filter(variantID %in% variants_select)

ASE_ASoC <- ASE_ASoC %>%
  mutate(ref_prop = refCount/totalCount)
data_summary <- ASE_ASoC %>%
  group_by(variantID, celltype, measurement, condition) %>%
  summarise(
    mean_prop = mean(ref_prop, na.rm = TRUE),
    n = n(),
    low = min(ref_prop, na.rm = TRUE),  
    high = max(ref_prop, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(n >= 1)
data_summary_sub <- data_summary %>%
  filter(variantID %in% c('5_40687361', '5_40687957', '5_40692838'))

data_summary_sub$celltype <- gsub("_", " ", data_summary_sub$celltype)
data_summary_sub$celltype <- factor(data_summary_sub$celltype, levels = c("CD4+ Naive Resting", "CD4+ Naive Activated", 
                                                                          "CD4+ Regulatory Resting", "CD4+ Regulatory Activated",
                                                                          "CD4+ Memory Resting Th1", "CD4+ Memory Activated Th1",
                                                                          "CD4+ Memory Resting Th17", "CD4+ Memory Activated Th17",
                                                                          "CD4+ Memory Resting Tfh", "CD4+ Memory Activated Tfh",
                                                                          "CD4+ Memory Resting Other", "CD4+ Memory Activated Other",
                                                                          "CD8+ Naive Resting", "CD8+ Naive Activated", 
                                                                          "CD8+ Regulatory", "CD8+ Memory Resting", 
                                                                          "CD8+ Memory Activated", "MAITs Resting", 
                                                                          "MAITs Activated", "Gamma Delta"))
ggplot(data_summary_sub, aes(x = mean_prop, 
                             y = celltype, color = measurement)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(xmin = low, xmax = high), 
                width = 0.3, 
                position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  xlim(0, 1) +
  scale_color_manual(values = c("ASE" = "#4DBBD5B2", "ASoC" = "#DC0000B2")) +
  facet_grid(variantID ~ condition, scales = "free_x") +
  labs(x = NULL, y = NULL, title = "Forest Plot of Conditions by Cell Type") +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 9),
        legend.position = "bottom")

#####
# ASE/ASoC overlaps
ASE_Results <- fread(paste0("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASE_Count/ASE_counts.txt")) %>%
  mutate(measurement = "ASE") %>% filter(p_value_adj < 0.05)
ASoC_Results <- fread(paste0("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASoC_Count/ASoC_counts.txt")) %>%
  mutate(measurement = "ASoC") %>% filter(p_value_adj < 0.05)
ASE_ASoC <- rbind(ASE_Results, ASoC_Results)
intersc <- intersect(ASE_Results$variantID, ASoC_Results$variantID)
ASE_ASoC <- ASE_ASoC %>% filter(variantID %in% intersc)
a <- table(ASE_ASoC$variantID) %>% 
  data.frame() %>%
  dplyr::arrange(desc(Freq))
