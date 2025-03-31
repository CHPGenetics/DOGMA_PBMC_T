rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(liftOver)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)

# Parameters
DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASoC_Count/"
SUMM_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Data/GWAS_Summary_Raw/"
REF_PATH <- "/ix1/wchen/Shiyue/References/LiftOver/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/"

# Load data
ASoC_Results <- fread(paste0(DATA_PATH, "ASoC_counts.txt"))
ASoC_Results <- ASoC_Results %>%
  mutate(sample_id = paste0(condition, "_", sample))
ASoC_Results$variantID <- paste0("chr", ASoC_Results$variantID)

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

# IBD GWAS
IBD <- fread(paste0(SUMM_PATH, "NG_IBD_3760/ibd_build37_59957_20161107.txt"))
IBD <- IBD %>%
  separate(col = MarkerName, into = c("CHR", "BP"), sep = ":") %>%
  mutate(BP = gsub("_.*", "", BP))
IBD <- data.frame(seqnames = as.numeric(IBD$CHR),
                  start = as.numeric(IBD$BP),
                  end = as.numeric(IBD$BP),
                  P = as.numeric(IBD$P.value),
                  A1 = toupper(IBD$Allele1),
                  A2 = toupper(IBD$Allele2))
## convert to hg38
## Liftover
cur <- makeGRangesFromDataFrame(IBD, keep.extra.columns = T)
ch = import.chain(paste0(REF_PATH, 'hg19ToHg38.over.chain'))
ch

str(ch[[1]])
seqlevelsStyle(cur) = "UCSC"  # necessary
cur38 = liftOver(cur, ch)
class(cur38)

cur38 = unlist(cur38)
genome(cur38) = "hg38"
cur38

IBD <- data.frame(CHR = seqnames(cur38),
                  BP = start(cur38),
                  P = cur38$P,
                  A1 = cur38$A1,
                  A2 = cur38$A2)

IBD <- IBD %>% 
  #  mutate(CHR = gsub("chr", "", CHR)) %>%
  mutate(variantID = paste0(CHR, "_", BP))
IBD <- IBD %>% filter(variantID %in% ASoC_Results$variantID) %>%
  mutate(condition = ifelse(variantID %in% PGE2_unique_snp, "Act_IL1B_IL23_PGE2-specific", 
                            ifelse(variantID %in% Control_unique_snp, "Act_IL1B_IL23-specific", 
                                   ifelse(variantID %in% TGFB_unique_snp, "Act_IL1B_IL23_TGFB-specific", 
                                          ifelse(variantID %in% TWO_unique_snp, 
                                                 "Act_IL1B_IL23_PGE2_TGFB-specific", "Shared")))))

# QQ plot
IBD <- IBD %>%
  arrange(P) %>%
  group_by(condition) %>%
  mutate(expected = -log10((rank(P, ties.method = "first") - 0.5) / length(P))) %>%
  ungroup() %>%
  mutate(observed = -log10(P)) %>%
  arrange(condition)

test_results <- lapply(unique(IBD$condition), function(cond) {
  IBD_sub <- IBD %>% filter(condition == cond)
  test <- wilcox.test(IBD_sub$observed, IBD_sub$expected)
  paste0(cond, ": Wilcoxon P = ", format(test$p.value, digits = 3))
})

test_labels <- paste(test_results, collapse = "\n")

pdf(paste0(WORK_PATH, "IBD_ASoC_overlap.pdf"))
ggplot() + 
  geom_point(data = IBD %>% filter(condition == "Shared"), 
             mapping = aes(x = expected, y = observed, color = condition), 
             alpha = 0.5) +
  geom_point(data = IBD %>% filter(condition != "Shared"), 
             mapping = aes(x = expected, y = observed, color = condition), 
             alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
  scale_color_manual(values = c("Shared" = "#B0B0B0FF", 
                                "Act_IL1B_IL23_PGE2-specific" = "#CD534CFF", 
                                "Act_IL1B_IL23_TGFB-specific" = "#0073C2FF", 
                                "Act_IL1B_IL23_PGE2_TGFB-specific" = "#EFC000FF", 
                                "Act_IL1B_IL23-specific" = "#00A087FF")) +
  labs(title = "IBD GWAS QQ Plot (ASoC signals)", x = "Expected -log10(p)", y = "Observed -log10(p)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  annotate("text", x = 0, y = 39.5, label = test_labels, 
           hjust = 0, vjust = 1, size = 4)
dev.off()


