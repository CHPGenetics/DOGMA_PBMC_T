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
DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASoC_Count/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/"

# Load data
ASoC_Results <- fread(paste0(DATA_PATH, "ASoC_counts.txt"))
ASoC_Results <- ASoC_Results %>%
  filter(p_value_adj < 0.05) %>%
  mutate(sample_id = paste0(condition, "_", sample))

# Summary
PGE2_ASoC <- ASoC_Results %>% filter(condition == "Act_IL1B_IL23_PGE2")
TGFB_ASoC <- ASoC_Results %>% filter(condition == "Act_IL1B_IL23_PGE2_TGFB")
TWO_ASoC <- ASoC_Results %>% filter(condition == "Act_IL1B_IL23_TGFB")
Control_ASoC <- ASoC_Results %>% filter(condition == "Act_IL1B_IL23")

# Proportion
PGE2_unique_snp <- setdiff(PGE2_ASoC$variantID, Control_ASoC$variantID) %>%
  setdiff(TGFB_ASoC$variantID) %>% setdiff(TWO_ASoC$variantID)
TGFB_unique_snp <- setdiff(TGFB_ASoC$variantID, Control_ASoC$variantID) %>%
  setdiff(PGE2_ASoC$variantID) %>% setdiff(TWO_ASoC$variantID)
Control_unique_snp <- setdiff(Control_ASoC$variantID, PGE2_ASoC$variantID) %>%
  setdiff(TGFB_ASoC$variantID) %>% setdiff(TWO_ASoC$variantID)
TWO_unique_snp <- setdiff(TWO_ASoC$variantID, PGE2_ASoC$variantID) %>%
  setdiff(TGFB_ASoC$variantID) %>% setdiff(Control_ASoC$variantID)
Shared_snp <- intersect(unique(PGE2_ASoC$variantID), unique(Control_ASoC$variantID)) %>%
  intersect(unique(TGFB_ASoC$variantID)) %>% intersect(unique(TWO_ASoC$variantID))


## Top ASoC cites
top_ASoC <- table(ASoC_Results$variantID) %>% 
  data.frame() %>%
  dplyr::arrange(desc(Freq))

## Top ASoC cites in Control
top_ASoC_Control <- ASoC_Results %>%
  filter(variantID %in% Control_unique_snp)
top_ASoC_Control <- table(top_ASoC_Control$variantID) %>% 
  data.frame() %>%
  dplyr::arrange(desc(Freq))

## Top ASoC cites in PGE2
top_ASoC_PGE2 <- ASoC_Results %>%
  filter(variantID %in% PGE2_unique_snp)
top_ASoC_PGE2 <- table(top_ASoC_PGE2$variantID) %>% 
  data.frame() %>%
  dplyr::arrange(desc(Freq))

## Top Shared ASoC cites
top_ASoC_shared <- ASoC_Results %>%
  filter(variantID %in% Shared_snp)
top_ASoC_shared <- table(top_ASoC_shared$variantID) %>% 
  data.frame() %>%
  dplyr::arrange(desc(Freq))

# Forest plot
ASoC_Results <- ASoC_Results %>%
  mutate(ref_prop = refCount/totalCount)

## Shared
ASoC_Results_sub <- ASoC_Results %>%
  filter(variantID %in% Shared_snp) %>%
  filter(variantID %in% top_ASoC$Var1[1:21]) %>%
  mutate(variantID = factor(variantID, levels = rev(top_ASoC$Var1))) %>%
  arrange(variantID)

pdf(paste0(WORK_PATH, "ASoC_Ref_proportion_shared.pdf"), width = 6, height = 5)
ggplot(ASoC_Results_sub, aes(y = variantID, x = ref_prop, color = condition)) +
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
ASoC_Results_sub <- ASoC_Results %>%
  filter(variantID %in% PGE2_unique_snp) %>%
  filter(variantID %in% top_ASoC_PGE2$Var1[1:10]) %>%
  mutate(variantID = factor(variantID, levels = rev(top_ASoC_PGE2$Var1))) %>%
  arrange(variantID)

pdf(paste0(WORK_PATH, "ASoC_Ref_proportion_PGE2.pdf"), width = 6, height = 5)
ggplot(ASoC_Results_sub, aes(y = variantID, x = ref_prop, color = condition)) +
  geom_point() +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  scale_color_manual(values = c("#DC0000B2"),
                     labels = c("Act_IL1B_IL23_PGE2")) +
  xlim(0, 1) +
  theme_bw() +
  labs(x = "Reference Allele Mapping Proportion", y = "SNPs")
dev.off()

## Unique (Control)
ASoC_Results_sub <- ASoC_Results %>%
  filter(variantID %in% Control_unique_snp) %>%
  filter(variantID %in% top_ASoC_Control$Var1[1:10]) %>%
  mutate(variantID = factor(variantID, levels = rev(top_ASoC_Control$Var1))) %>%
  arrange(variantID)

pdf(paste0(WORK_PATH, "ASoC_Ref_proportion_Control.pdf"), width = 6, height = 5)
ggplot(ASoC_Results_sub, aes(y = variantID, x = ref_prop, color = condition)) +
  geom_point() +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  scale_color_manual(values = c("#4DBBD5B2"),
                     labels = c("Act_IL1B_IL23")) +
  xlim(0, 1) +
  theme_bw() +
  labs(x = "Reference Allele Mapping Proportion", y = "SNPs")
dev.off()

# Upset plot
input <- c(Control = length(unique(ASoC_Results$variantID[which(ASoC_Results$condition == "Act_IL1B_IL23")])),
           PGE2 = length(unique(ASoC_Results$variantID[which(ASoC_Results$condition == "Act_IL1B_IL23_")])),
           TGFB = length(unique(ASoC_Results$variantID[which(ASoC_Results$condition == "Act_IL1B_IL23_TGFB")])),
           PGE2_TGFB = length(unique(ASoC_Results$variantID[which(ASoC_Results$condition == "Act_IL1B_IL23_PGE2_TGFB")])),
           "Control&PGE2" = length(intersect(unique(PGE2_ASoC$variantID), unique(Control_ASoC$variantID))),
           "Control&TGFB" = length(intersect(unique(TGFB_ASoC$variantID), unique(Control_ASoC$variantID))),
           "Control&PGE2_TGFB" = length(intersect(unique(TWO_ASoC$variantID), unique(Control_ASoC$variantID))),
           "PGE2&TGFB" = length(intersect(unique(PGE2_ASoC$variantID), unique(TGFB_ASoC$variantID))),
           "PGE2&PGE2_TGFB" = length(intersect(unique(PGE2_ASoC$variantID), unique(TWO_ASoC$variantID))),
           "TGFB&PGE2_TGFB" = length(intersect(unique(TGFB_ASoC$variantID), unique(TWO_ASoC$variantID))),
           "Control&PGE2&TGFB&PGE2_TGFB" = length(Shared_snp))

pdf(paste0(WORK_PATH, "ASoC_Upset_plot.pdf"), width = 12, height = 6)
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
x <- list(Act_IL1B_IL23_PGE2 = unique(ASoC_Results$variantID[which(ASoC_Results$condition == "Act_IL1B_IL23_PGE2")]),
          Act_IL1B_IL23_PGE2_TGFB = unique(ASoC_Results$variantID[which(ASoC_Results$condition == "Act_IL1B_IL23_PGE2_TGFB")]),
          Act_IL1B_IL23_TGFB = unique(ASoC_Results$variantID[which(ASoC_Results$condition == "Act_IL1B_IL23_TGFB")]),
          Act_IL1B_IL23 = unique(ASoC_Results$variantID[which(ASoC_Results$condition == "Act_IL1B_IL23")]))

pdf(paste0(WORK_PATH, "ASoC_Venn_plot.pdf"), width = 11, height = 6.8)
ggvenn(x, 
       fill_color = c("#CD534CFF", "#0073C2FF", "#EFC000FF", "#868686FF"),
       stroke_size = 0.5, set_name_size = 3) # Reduce set label size
dev.off()

