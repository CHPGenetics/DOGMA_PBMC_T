library(data.table)
library(stringr)
library(liftOver)
library(SNPlocs.Hsapiens.dbSNP141.GRCh38)
library(rtracklayer)
library(Signac)
library(cowplot)
library(ggsci)
library(ggpubr)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

read.dir <- function(path, pattern){
  files <- list.files(path, pattern)
  tables <- list()
  for (i in files){
    tmp <- read.table(paste(path, i, sep = '/'), sep = '\t')
    colnames(tmp) <- c('name', 'p_value')
    tmp$SNP <- gsub(pattern, '', i)
    tables[[gsub(pattern, '', i)]] <- tmp
  }
  return(tables)
}

table <- read.dir('../output/SNP/CHEERS/AllPeaks/Merged_analysis/', '_disease_enrichment_pValues.txt')

data_plot_bp <- table[[1]]
for(i in 2:length(table)){
  data_plot_bp <- rbind(data_plot_bp, table[[i]])
}

data_plot <- data_plot_bp[grep('IBD|UC|CD', data_plot_bp$SNP, invert = T),]
data_plot <- data_plot[grep('Immune|Non-immune', data_plot$SNP),]

pdf('../plots/SNP/CHEERS/enrich_immune.pdf', width = 4, height = 2)
ggplot(data_plot, aes(y = -log10(p_value), x = SNP, fill = SNP))+
  geom_violin()+
  geom_point(color = 'black', shape = 21)+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'darkgrey')+
  scale_fill_manual(values = c('firebrick', 'dodgerblue3'))+
  coord_flip()+
  theme_cowplot()+
  theme(legend.position = 'none')+
  labs(x = 'SNPs', y = '-log10(P-value)')
dev.off()

pdf('../plots/SNP/CHEERS/enrich_immune_paired.pdf', width = 4, height = 2)
ggplot(data_plot, aes(y = -log10(p_value), x = SNP, col = SNP))+
  geom_violin(fill = 'white',)+
  geom_point(fill = 'white', shape = 21)+
  geom_line(aes(group = name), col = 'darkgrey')+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'darkgrey')+
  scale_color_manual(values = c('firebrick', 'dodgerblue3'))+
  coord_flip()+
  theme_cowplot()+
  theme(legend.position = 'none')+
  labs(x = 'SNPs', y = '-log10(P-value)')
dev.off()

data_plot <- data_plot_bp[grep('IBD|UC|CD', data_plot_bp$SNP, invert = T),]
data_plot <- data_plot[grep('Immune|Non-immune', data_plot$SNP, invert = T),]
data_plot$Group <- 'Immune'
data_plot$Group[data_plot$SNP %in% c('Platelet_counts', 'Red_blood_cell_traits',
                                     'Urate_levels', 'Triglycerides', 'C_reactive_protein',
                                     'HDL_cholesterol', 'LDL_cholesterol',
                                     'Liver_enzyme_levels_gamma_glutamyl_transferase', 'Creatinine_levels')] <- 'Non-immune 1'
data_plot$Group[data_plot$SNP %in% c('Migraine', 'Renal_function_related_traits_BUN', 'Chronic_kidney_disease',
                                     'Type_2_diabetes', 'Bone_mineral_density', 'Fasting_glucose_related_traits',
                                     'Alzheimers_combined', 'Restless_legs_syndrome', 'Progressive_supranuclear_palsy')] <- 'Non-immune 2'
data_plot$log10_p_value <- -log10(data_plot$p_value)
data_plot$log10_p_value[data_plot$p_value >= 0.05] <- NA

data_plot1 <- data_plot[data_plot$Group == 'Immune',]
data_plot1$SNP <- factor(data_plot1$SNP, levels = names(table(data_plot1$SNP)))
data_plot2 <- data_plot[data_plot$Group == 'Non-immune 1',]
data_plot2$SNP <- factor(data_plot2$SNP, levels = names(table(data_plot2$SNP)))
data_plot3 <- data_plot[data_plot$Group == 'Non-immune 2',]
data_plot3$SNP <- factor(data_plot3$SNP, levels = names(table(data_plot3$SNP)))

p1 <- ggplot(data_plot1, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,5), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)')+ 
  coord_fixed(expand = T)
p2 <- ggplot(data_plot2, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,5), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(P-value)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
p3 <- ggplot(data_plot3, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,5), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/CHEERS/enrich_disease.pdf', width = 12, height = 12)
p3/p2/p1
dev.off()

p1 <- ggplot(data_plot1, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,5), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(P-value)')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/CHEERS/enrich_disease_immune.pdf', width = 8, height = 8)
p1
dev.off()

############
data_plot <- data_plot_bp[grep('IBD|UC|CD', data_plot_bp$SNP),]
data_plot <- data_plot[grep('IBD', data_plot$SNP),]
data_plot$Group <- 'All variants'
data_plot$Group[grep('Lead', data_plot$SNP)] <- 'Leading variants'
data_plot$Group[grep('Causal', data_plot$SNP)] <- 'Causal variants with > 95% certainty'

pdf('../plots/SNP/CHEERS/enrich_immune_IBD.pdf', width = 6, height = 3)
ggplot(data_plot, aes(y = -log10(p_value), x = Group, fill = Group))+
  geom_violin()+
  geom_point(color = 'black', shape = 21)+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'darkgrey')+
  scale_fill_npg()+
  coord_flip()+
  theme_cowplot()+
  theme(legend.position = 'none')+
  labs(x = 'SNPs', y = '-log10(P-value)')
dev.off()

pdf('../plots/SNP/CHEERS/enrich_immune_paired.pdf', width = 6, height = 3)
ggplot(data_plot, aes(y = -log10(p_value), x = Group, col = Group))+
  geom_violin(fill = 'white',)+
  geom_point(fill = 'white', shape = 21)+
  geom_line(aes(group = name), col = 'darkgrey')+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'darkgrey')+
  scale_color_npg()+
  coord_flip()+
  theme_cowplot()+
  theme(legend.position = 'none')+
  labs(x = 'SNPs', y = '-log10(P-value)')
dev.off()

data_plot <- data_plot_bp[grep('IBD|UC|CD', data_plot_bp$SNP),]
data_plot$snp <- 'All variants'
data_plot$snp[grep('Lead', data_plot$SNP)] <- 'Leading variants'
data_plot$snp[grep('Causal', data_plot$SNP)] <- 'Causal variants with > 95% certainty'
data_plot$snp <- factor(data_plot$snp, levels = c('Causal variants with > 95% certainty', 'Leading variants', 'All variants'))
data_plot$group <- 'IBD'
data_plot$group[grep('UC', data_plot$SNP)] <- 'UC'
data_plot$group[grep('CD', data_plot$SNP)] <- 'CD'
data_plot$snp_group <- paste0(data_plot$snp, ' (', data_plot$group, ')')

data_plot$log10_p_value <- -log10(data_plot$p_value)
data_plot$log10_p_value[data_plot$p_value >= 0.05] <- NA

data_plot1 <- data_plot[data_plot$group == 'IBD',]
data_plot1$snp_group <- factor(data_plot1$snp_group, levels = c('Causal variants with > 95% certainty (IBD)', 'Leading variants (IBD)', 'All variants (IBD)'))
data_plot2 <- data_plot[data_plot$group == 'UC',]
data_plot2$snp_group <- factor(data_plot2$snp_group, levels = c('Causal variants with > 95% certainty (UC)', 'Leading variants (UC)', 'All variants (UC)'))
data_plot3 <- data_plot[data_plot$group == 'CD',]
data_plot3$snp_group <- factor(data_plot3$snp_group, levels = c('Causal variants with > 95% certainty (CD)', 'Leading variants (CD)', 'All variants (CD)'))

p1 <- ggplot(data_plot1, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,5), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'bottom')+
  labs(x = '', y = '', fill = '-log10(P-value)')+ 
  coord_fixed(expand = T)
p2 <- ggplot(data_plot2, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,5), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
p3 <- ggplot(data_plot3, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,5), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/CHEERS/enrich_disease_IBD.pdf', width = 10, height = 8)
p3/p2/p1
dev.off()

###########
data_plot$snp <- gsub('Leading', 'Lead', data_plot$snp)
data_plot$snp_group <- gsub('Leading', 'Lead', data_plot$snp_group)
data_plot <- data_plot[data_plot$snp != 'Causal variants with > 95% certainty',]
data_plot1 <- data_plot[data_plot$group == 'IBD',]
data_plot1$snp_group <- factor(data_plot1$snp_group, levels = c('Lead variants (IBD)', 'All variants (IBD)'))
data_plot2 <- data_plot[data_plot$group == 'UC',]
data_plot2$snp_group <- factor(data_plot2$snp_group, levels = c('Lead variants (UC)', 'All variants (UC)'))
data_plot3 <- data_plot[data_plot$group == 'CD',]
data_plot3$snp_group <- factor(data_plot3$snp_group, levels = c('Lead variants (CD)', 'All variants (CD)'))

p1 <- ggplot(data_plot1, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,4), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)')+ 
  coord_fixed(expand = T)
p2 <- ggplot(data_plot2, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,4), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(P-value)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
p3 <- ggplot(data_plot3, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,4), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/CHEERS/enrich_disease_IBD_updated.pdf', width = 10, height = 6.5)
p3/p2/p1
dev.off()
######################
table <- read.dir('../output/SNP/CHEERS/AllPeaks/Individual_analysis/', '_disease_enrichment_pValues.txt')

data_plot_bp <- table[[1]]
for(i in 2:length(table)){
  data_plot_bp <- rbind(data_plot_bp, table[[i]])
}
data_plot_bp$celltype <- paste(str_split_fixed(data_plot_bp$name, '_', 3)[,1], str_split_fixed(data_plot_bp$name, '_', 3)[,2], sep = '_')
data_plot_bp$sample <- str_split_fixed(data_plot_bp$name, '_', 3)[,3]
data_plot_bp$sample <- as.factor(data_plot_bp$sample)

# data_plot <- data_plot_bp[grep('IBD|UC|CD', data_plot_bp$SNP, invert = T),]
# data_plot <- data_plot[grep('Immune|Non-immune', data_plot$SNP),]
# 
# pdf('../plots/SNP/CHEERS/enrich_immune_individual.pdf', width = 8, height = 8)
# ggplot(data_plot, aes(y = -log10(p_value), x = SNP, fill = SNP))+
#   geom_boxplot()+
#   geom_point(color = 'black', shape = 21)+
#   geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'darkgrey')+
#   scale_fill_manual(values = c('firebrick', 'dodgerblue3'))+
#   coord_flip()+
#   facet_wrap(~celltype, ncol = 1, strip.position="left")+
#   theme_cowplot()+
#   theme(strip.background = element_rect(fill = 'white'), strip.text.y.left = element_text(angle = 0))+
#   labs(x = '', y = '-log10(P-value)', fill = 'SNPs')+
#   rremove('y.text')
# dev.off()
# 
# pdf('../plots/SNP/CHEERS/enrich_immune_paired_individual.pdf', width = 8, height = 8)
# ggplot(data_plot, aes(y = -log10(p_value), x = SNP, col = SNP))+
#   geom_boxplot(fill = 'white')+
#   geom_point(fill = 'white', shape = 21)+
#   geom_line(aes(group = name), col = 'darkgrey')+
#   geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'darkgrey')+
#   scale_color_manual(values = c('firebrick', 'dodgerblue3'))+
#   coord_flip()+
#   facet_wrap(~celltype, ncol = 1, strip.position="left")+
#   theme_cowplot()+
#   theme(strip.background = element_rect(fill = 'white'), strip.text.y.left = element_text(angle = 0))+
#   labs(x = '', y = '-log10(P-value)', color = 'SNPs')+
#   rremove('y.text')
# dev.off()

data_plot <- data_plot_bp[grep('IBD|UC|CD', data_plot_bp$SNP, invert = T),]
data_plot <- data_plot[grep('Immune|Non-immune', data_plot$SNP, invert = T),]
data_plot$Group <- 'Immune'
data_plot$Group[data_plot$SNP %in% c('Platelet_counts', 'Red_blood_cell_traits',
                                     'Urate_levels', 'Triglycerides', 'C_reactive_protein',
                                     'HDL_cholesterol', 'LDL_cholesterol',
                                     'Liver_enzyme_levels_gamma_glutamyl_transferase', 'Creatinine_levels')] <- 'Non-immune 1'
data_plot$Group[data_plot$SNP %in% c('Migraine', 'Renal_function_related_traits_BUN', 'Chronic_kidney_disease',
                                     'Type_2_diabetes', 'Bone_mineral_density', 'Fasting_glucose_related_traits',
                                     'Alzheimers_combined', 'Restless_legs_syndrome', 'Progressive_supranuclear_palsy')] <- 'Non-immune 2'
# data_plot$log10_p_value <- -log10(data_plot$p_value)
# data_plot$log10_p_value[data_plot$p_value >= 0.05] <- NA

data_plot1 <- data_plot[data_plot$Group == 'Immune',]
data_plot1$SNP <- factor(data_plot1$SNP, levels = names(table(data_plot1$SNP)))
data_plot2 <- data_plot[data_plot$Group == 'Non-immune 1',]
data_plot2$SNP <- factor(data_plot2$SNP, levels = names(table(data_plot2$SNP)))
data_plot3 <- data_plot[data_plot$Group == 'Non-immune 2',]
data_plot3$SNP <- factor(data_plot3$SNP, levels = names(table(data_plot3$SNP)))

df_tmp = data_plot1[data_plot1$SNP == 'Allergy',]
pdf('../plots/SNP/CHEERS/enrich_disease_individual_legend.pdf', width = 8, height = 8)
ggplot(df_tmp, aes(y = -log10(p_value), x = celltype, fill = celltype))+
  geom_boxplot()+
  geom_point(color = 'black', shape = 21)+
  scale_fill_discrete(breaks=names(table(df_tmp$celltype))[20:1])+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'darkgrey')+
  # scale_fill_manual(values = c('firebrick', 'dodgerblue3'))+
  coord_flip()+
  theme_cowplot()+
  theme(legend.position = 'right')+
  labs(x = 'Cell types', y = '-log10(P-value)', fill = '')+
  rremove('y.text')+
  guides(fill = guide_legend(ncol = 1))
dev.off()

for(i in 1:21){
  name <- names(table(data_plot1$SNP))[i]
  assign(paste0('p',i), ggplot(data_plot1[data_plot1$SNP == name,], aes(y = -log10(p_value), x = celltype, fill = celltype))+
           geom_boxplot()+
           geom_point(color = 'black', shape = 21)+
           geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'firebrick')+
           # scale_fill_manual(values = c('firebrick', 'dodgerblue3'))+
           coord_flip()+
           theme_cowplot()+
           ggtitle(name)+
           theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
           labs(x = 'Cell types', y = '-log10(P-value)', fill = 'SNPs')+
           rremove('y.text'))
}

pdf('../plots/SNP/CHEERS/enrich_disease_individual_immune.pdf', width = 28, height = 12)
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21, nrow = 3)
dev.off()

pdf('../plots/SNP/CHEERS/enrich_disease_individual_immune_1.pdf', width = 20, height = 12)
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15, nrow = 3)
dev.off()

pdf('../plots/SNP/CHEERS/enrich_disease_individual_immune_2.pdf', width = 8, height = 12)
wrap_plots(p16,p17,p18,p19,p20,p21, nrow = 3)
dev.off()

for(i in 1:9){
  name <- names(table(data_plot2$SNP))[i]
  assign(paste0('p',i), ggplot(data_plot2[data_plot2$SNP == name,], aes(y = -log10(p_value), x = celltype, fill = celltype))+
           geom_boxplot()+
           geom_point(color = 'black', shape = 21)+
           geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'firebrick')+
           # scale_fill_manual(values = c('firebrick', 'dodgerblue3'))+
           coord_flip()+
           theme_cowplot()+
           ggtitle(name)+
           theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
           labs(x = 'Cell types', y = '-log10(P-value)', fill = 'SNPs')+
           rremove('y.text'))
}

pdf('../plots/SNP/CHEERS/enrich_disease_individual_nonimmune1.pdf', width = 12, height = 12)
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8,p9, nrow = 3)
dev.off()

for(i in 1:9){
  name <- names(table(data_plot3$SNP))[i]
  assign(paste0('p',i), ggplot(data_plot3[data_plot3$SNP == name,], aes(y = -log10(p_value), x = celltype, fill = celltype))+
           geom_boxplot()+
           geom_point(color = 'black', shape = 21)+
           geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'firebrick')+
           # scale_fill_manual(values = c('firebrick', 'dodgerblue3'))+
           coord_flip()+
           theme_cowplot()+
           ggtitle(name)+
           theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
           labs(x = 'Cell types', y = '-log10(P-value)', fill = 'SNPs')+
           rremove('y.text'))
}

pdf('../plots/SNP/CHEERS/enrich_disease_individual_nonimmune2.pdf', width = 12, height = 12)
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8,p9, nrow = 3)
dev.off()

############
# data_plot <- data_plot_bp[grep('IBD|UC|CD', data_plot_bp$SNP),]
# data_plot <- data_plot[grep('IBD', data_plot$SNP),]
# data_plot$Group <- 'All variants'
# data_plot$Group[grep('Lead', data_plot$SNP)] <- 'Leading variants'
# data_plot$Group[grep('Causal', data_plot$SNP)] <- 'Causal variants with > 95% certainty'
# 
# pdf('../plots/SNP/CHEERS/enrich_immune_IBD.pdf', width = 4, height = 2)
# ggplot(data_plot, aes(y = -log10(p_value), x = Group, fill = Group))+
#   geom_violin()+
#   geom_point(color = 'black', shape = 21)+
#   geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'darkgrey')+
#   scale_fill_npg()+
#   coord_flip()+
#   theme_cowplot()+
#   theme(legend.position = 'none')+
#   labs(x = 'SNPs', y = '-log10(P-value)')
# dev.off()
# 
# pdf('../plots/SNP/CHEERS/enrich_immune_paired_IBD.pdf', width = 4, height = 2)
# ggplot(data_plot, aes(y = -log10(p_value), x = Group, col = Group))+
#   geom_violin(fill = 'white',)+
#   geom_point(fill = 'white', shape = 21)+
#   geom_line(aes(group = name), col = 'darkgrey')+
#   geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'darkgrey')+
#   scale_color_npg()+
#   coord_flip()+
#   theme_cowplot()+
#   theme(legend.position = 'none')+
#   labs(x = 'SNPs', y = '-log10(P-value)')
# dev.off()

data_plot <- data_plot_bp[grep('IBD|UC|CD', data_plot_bp$SNP),]
data_plot$snp <- 'All variants'
data_plot$snp[grep('Lead', data_plot$SNP)] <- 'Leading variants'
data_plot$snp[grep('Causal', data_plot$SNP)] <- 'Causal variants with > 95% certainty'
data_plot$snp <- factor(data_plot$snp, levels = c('Causal variants with > 95% certainty', 'Leading variants', 'All variants'))
data_plot$group <- 'IBD'
data_plot$group[grep('UC', data_plot$SNP)] <- 'UC'
data_plot$group[grep('CD', data_plot$SNP)] <- 'CD'
data_plot$snp_group <- paste0(data_plot$snp, ' (', data_plot$group, ')')

# data_plot$log10_p_value <- -log10(data_plot$p_value)
# data_plot$log10_p_value[data_plot$p_value >= 0.1] <- NA

data_plot1 <- data_plot[data_plot$group == 'IBD',]
data_plot1$snp_group <- factor(data_plot1$snp_group, levels = c('Causal variants with > 95% certainty (IBD)', 'Leading variants (IBD)', 'All variants (IBD)'))
data_plot2 <- data_plot[data_plot$group == 'UC',]
data_plot2$snp_group <- factor(data_plot2$snp_group, levels = c('Causal variants with > 95% certainty (UC)', 'Leading variants (UC)', 'All variants (UC)'))
data_plot3 <- data_plot[data_plot$group == 'CD',]
data_plot3$snp_group <- factor(data_plot3$snp_group, levels = c('Causal variants with > 95% certainty (CD)', 'Leading variants (CD)', 'All variants (CD)'))
data_plot$snp_group <- factor(data_plot$snp_group, levels = c('Causal variants with > 95% certainty (IBD)', 'Leading variants (IBD)', 'All variants (IBD)',
                                                              'Causal variants with > 95% certainty (UC)', 'Leading variants (UC)', 'All variants (UC)',
                                                              'Causal variants with > 95% certainty (CD)', 'Leading variants (CD)', 'All variants (CD)'))

for(i in 1:9){
  name <- names(table(data_plot$snp_group))[i]
  assign(paste0('p',i), ggplot(data_plot[data_plot$snp_group == name,], aes(y = -log10(p_value), x = celltype, fill = celltype))+
           geom_boxplot()+
           geom_point(color = 'black', shape = 21)+
           geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'firebrick')+
           # scale_fill_manual(values = c('firebrick', 'dodgerblue3'))+
           coord_flip()+
           theme_cowplot()+
           ggtitle(name)+
           theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
           labs(x = 'Cell types', y = '-log10(P-value)', fill = 'SNPs')+
           rremove('y.text'))
}

pdf('../plots/SNP/CHEERS/enrich_disease_IBD_individual.pdf', width = 12, height = 12)
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8,p9, nrow = 3)
dev.off()

###########
data_plot$snp <- gsub('Leading', 'Lead', data_plot$snp)
data_plot$snp_group <- gsub('Leading', 'Lead', data_plot$snp_group)
data_plot <- data_plot[data_plot$snp != 'Causal variants with > 95% certainty',]
data_plot$snp_group <- factor(as.character(data_plot$snp_group), levels = c('Lead variants (IBD)', 'All variants (IBD)',
                                                              'Lead variants (UC)', 'All variants (UC)',
                                                              'Lead variants (CD)', 'All variants (CD)'))



for(i in 1:6){
  name <- names(table(data_plot$snp_group))[i]
  assign(paste0('p',i), ggplot(data_plot[data_plot$snp_group == name,], aes(y = -log10(p_value), x = celltype, fill = celltype))+
           geom_boxplot()+
           geom_point(color = 'black', shape = 21)+
           geom_hline(yintercept = -log10(0.05), linetype = 'dashed', col = 'firebrick')+
           # scale_fill_manual(values = c('firebrick', 'dodgerblue3'))+
           coord_flip()+
           theme_cowplot()+
           ggtitle(name)+
           theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
           labs(x = 'Cell types', y = '-log10(P-value)', fill = 'SNPs')+
           rremove('y.text'))
}

pdf('../plots/SNP/CHEERS/enrich_disease_IBD_individual_updated.pdf', width = 8, height = 12)
wrap_plots(p1,p2,p3,p4,p5,p6, nrow = 3)
dev.off()

 