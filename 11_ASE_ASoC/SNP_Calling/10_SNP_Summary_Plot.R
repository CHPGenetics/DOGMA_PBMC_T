rm(list=ls())
gc()

# Load Packages
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggsci)

# Parameters
DATA_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/03_SNP/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/Result/04_SNP_Summary/"

summary <- fread(paste0(WORK_PATH, "SNP_Number_overall.txt"))
proportion <- fread(paste0(WORK_PATH, "SNP_Proportion_overall.txt"))
# SNP Number
## Bar plot
summary$MAF <- as.character(summary$MAF)
summary$Depth <- as.character(summary$Depth)
summary_subset <- summary %>%
  filter(KGP == "KGP")

p <- ggplot(summary_subset, aes(Depth, N, fill = MAF)) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = N),                
            position = position_dodge(0.9),         
            vjust = -0.5) +                        
  facet_wrap(~ Resource, scales = "free") +
  labs(fill = "MAF", x = "Depth", y = "Number of SNPs") + 
  scale_y_continuous(labels = scales::comma,
                     limits = c(0, max(summary_subset$N) + 10)) +
  theme_bw() +
  scale_fill_manual(values = c("#7a929e", "#666464", "#68ac57", "#a8c498"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), 
        strip.background = element_blank(),
        axis.line.y.left = element_line(size = 1, color = "#E0E0E0"), 
        axis.line.x = element_line(size = 1, color = "#E0E0E0"))
pdf(paste0(WORK_PATH, "SNP_Number_Bar_Plot_Depth_KGP_Overall.pdf"), width = 32, height = 5)
print(p)
dev.off()

summary_subset <- summary %>%
  filter(KGP %in% c("Geno 0.5", "Geno KGP"))

p <- ggplot(summary_subset, aes(Depth, N, fill = KGP)) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = N),                
            position = position_dodge(0.9),         
            vjust = -0.5) +                        
  facet_wrap(~ Resource, scales = "free") +
  labs(fill = "KGP", x = "Depth", y = "Number of SNPs") + 
  scale_y_continuous(labels = scales::comma,
                     limits = c(0, max(summary_subset$N) + 10)) +
  theme_bw() +
  scale_fill_manual(values = c("#7a929e", "#68ac57"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), 
        strip.background = element_blank(),
        axis.line.y.left = element_line(size = 1, color = "#E0E0E0"), 
        axis.line.x = element_line(size = 1, color = "#E0E0E0"))

pdf(paste0(WORK_PATH, "SNP_Number_Bar_Plot_Depth_Geno_Overall.pdf"), width = 20, height = 5)
print(p)
dev.off()

## Scatter plot
summary$Depth <- as.character(summary$Depth)
p <- ggplot(summary, aes(x = MAF, y = N, color = Depth)) +
  geom_point(size = 3, alpha = 0.7) +              
  facet_wrap(~ KGP + Resource, scales = "free") + 
  scale_y_continuous(labels = scales::comma) +
  scale_color_jama() +                                  
  theme_classic() +                                     
  labs(x = "Minor Allele Frequency (MAF)",
       y = "SNP Number",
       color = "Resource") +                           
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),  
        legend.position = "top")
pdf(paste0(WORK_PATH, "Scatter_Plot_overall.pdf"), width = 8, height = 6)
print(p)
dev.off()

# MAF distribution hist plot
ref <- fread("/ix1/wchen/Shiyue/References/1KGP/hg38/frq_all.txt")
hist_sum <- data.frame()
for(depth in c(10, 20, 30, 40)){
  
  ## ATAC
  bim <- fread(paste0(DATA_PATH, "Combine/Final/ATAC_", depth, ".bim"))
  ref_a <- ref %>%
    filter(ID %in% bim$V2) %>%
    mutate(Depth = depth) %>%
    mutate(Resource = "ATAC")
  
  ## RNA
  bim <- fread(paste0(DATA_PATH, "Combine/Final/RNA_", depth, ".bim"))
  ref_b <- ref %>%
    filter(ID %in% bim$V2) %>%
    mutate(Depth = depth) %>%
    mutate(Resource = "RNA")
  
  ## Merge
  bim <- fread(paste0(DATA_PATH, "Combine/Final/merged_", depth, ".bim"))
  ref_c <- ref %>%
    filter(ID %in% bim$V2) %>%
    mutate(Depth = depth) %>%
    mutate(Resource = "Merge")
  
  hist_sum <- rbind(ref_a, ref_b, ref_c)
  p <- ggplot(hist_sum, aes(x = MAF)) +
    geom_histogram(bins = 50,   
                   fill = "#7a929e",
                   color = "black") + 
    facet_wrap(~ Resource, scales = "free") + 
    theme_classic() +                   
    labs(title = paste0("Depth: ", depth),
         x = "Minor Allele Frequency (MAF)", 
         y = "Density") +                
    theme(text = element_text(size = 14))
  
  pdf(paste0(WORK_PATH, "Hist_MAF_", depth, ".pdf"), width = 12, height = 6)
  print(p)
  dev.off()
}


#######################################################################################################
# Proportion plot
proportion$Prop <- round(proportion$Prop*100, 2)
proportion$Depth <- as.character(proportion$Depth)
p <- ggplot(proportion, aes(Type, Prop, fill = Depth)) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Proportion),                
            position = position_dodge(0.9),         
            vjust = -0.5) +                        
  facet_wrap(~ Resource, scales = "free") +
  labs(fill = "Depth", x = "", y = "Number of SNPs") + 
  scale_y_continuous(labels = scales::comma,
                     limits = c(0, max(proportion$Prop) + 10)) +
  theme_bw() +
  scale_fill_manual(values = c("#7a929e", "#666464", "#68ac57", "#a8c498"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), 
        strip.background = element_blank(),
        axis.line.y.left = element_line(size = 1, color = "#E0E0E0"), 
        axis.line.x = element_line(size = 1, color = "#E0E0E0"))

pdf(paste0(WORK_PATH, "SNP_Accuracy_overall.pdf"), width = 18, height = 6)
print(p)
dev.off()
