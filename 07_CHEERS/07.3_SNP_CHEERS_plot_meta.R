library(data.table)
library(stringr)
library(liftOver)
library(SNPlocs.Hsapiens.dbSNP141.GRCh38)
library(rtracklayer)
library(Signac)
library(cowplot)
library(ggsci)
library(ggpubr)
library(viridis)

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

table <- read.dir('../output/SNP/CHEERS/AllPeaks/Individual_analysis/', '_disease_enrichment_pValues.txt')

data_plot_bp <- table[[1]]
for(i in 2:length(table)){
  data_plot_bp <- rbind(data_plot_bp, table[[i]])
}
data_plot_bp$celltype <- paste(str_split_fixed(data_plot_bp$name, '_', 3)[,1], str_split_fixed(data_plot_bp$name, '_', 3)[,2], sep = '_')
data_plot_bp$sample <- str_split_fixed(data_plot_bp$name, '_', 3)[,3]
data_plot_bp$sample <- as.factor(data_plot_bp$sample)

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


meta_sample = function(df, pheno){
    df0 = df[df$SNP == pheno,]

    meta_p = sapply(names(table(df0$celltype)), function(i){
        df00 = df0[df0$celltype == i,]
        meta_p = metap::sumlog(df00$p_value)
        return(meta_p$p)
    })
    
    table = data.frame(name = names(meta_p), p_value = meta_p, SNP = pheno, Group = names(table(df0$Group)))
    return(table)
}

meta_list = lapply(names(table(data_plot$SNP)), function(pheno){
    meta_sample(data_plot, pheno)
})

meta = Reduce(rbind, meta_list)

data_plot = meta

data_plot$log10_p_value <- -log10(data_plot$p_value)
hist(data_plot$log10_p_value)
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
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)')+ 
  coord_fixed(expand = T)
p2 <- ggplot(data_plot2, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(P-value)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
p3 <- ggplot(data_plot3, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/CHEERS/enrich_disease_individual_meta.pdf', width = 12, height = 12)
p3/p2/p1
dev.off()

p1 <- ggplot(data_plot1, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/CHEERS/enrich_disease_individual_meta_immune.pdf', width = 8, height = 8)
p1
dev.off()
################
meta_list_fdr = meta_list
for (i in 1:39){
    meta_list_fdr[[i]]$p_val_adj = p.adjust(meta_list_fdr[[i]]$p_value, method = 'fdr')
}

meta = Reduce(rbind, meta_list_fdr)

data_plot = meta

data_plot$log10_p_value <- -log10(data_plot$p_val_adj)
hist(data_plot$log10_p_value)
data_plot$log10_p_value[data_plot$p_val_adj >= 0.05] <- NA

data_plot1 <- data_plot[data_plot$Group == 'Immune',]
data_plot1$SNP <- factor(data_plot1$SNP, levels = names(table(data_plot1$SNP)))
data_plot2 <- data_plot[data_plot$Group == 'Non-immune 1',]
data_plot2$SNP <- factor(data_plot2$SNP, levels = names(table(data_plot2$SNP)))
data_plot3 <- data_plot[data_plot$Group == 'Non-immune 2',]
data_plot3$SNP <- factor(data_plot3$SNP, levels = names(table(data_plot3$SNP)))

p1 <- ggplot(data_plot1, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(FDR)')+ 
  coord_fixed(expand = T)
p2 <- ggplot(data_plot2, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(FDR)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
p3 <- ggplot(data_plot3, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(FDR)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/CHEERS/enrich_disease_individual_meta_fdr.pdf', width = 12, height = 12)
p3/p2/p1
dev.off()

p1 <- ggplot(data_plot1, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(FDR)')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/CHEERS/enrich_disease_individual_meta_fdr_immune.pdf', width = 8, height = 8)
p1
dev.off()

################
meta_list_fdr = meta_list
for (i in 1:39){
    meta_list_fdr[[i]]$p_val_adj = p.adjust(meta_list_fdr[[i]]$p_value, method = 'bonferroni')
}

meta = Reduce(rbind, meta_list_fdr)
fwrite(meta, file = '../plots/SNP/CHEERS/enrich_disease_individual_meta_bonferroni.txt', col.names = T, row.names = F)

data_plot = meta

data_plot$log10_p_value <- -log10(data_plot$p_val_adj)
hist(data_plot$log10_p_value)
data_plot$log10_p_value[data_plot$p_val_adj >= 0.05] <- NA

data_plot1 <- data_plot[data_plot$Group == 'Immune',]
data_plot1$SNP <- factor(data_plot1$SNP, levels = names(table(data_plot1$SNP)))
data_plot2 <- data_plot[data_plot$Group == 'Non-immune 1',]
data_plot2$SNP <- factor(data_plot2$SNP, levels = names(table(data_plot2$SNP)))
data_plot3 <- data_plot[data_plot$Group == 'Non-immune 2',]
data_plot3$SNP <- factor(data_plot3$SNP, levels = names(table(data_plot3$SNP)))

p1 <- ggplot(data_plot1, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(Bonferroni)')+ 
  coord_fixed(expand = T)
p2 <- ggplot(data_plot2, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(Bonferroni)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
p3 <- ggplot(data_plot3, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(Bonferroni)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/CHEERS/enrich_disease_individual_meta_bonferroni.pdf', width = 12, height = 12)
p3/p2/p1
dev.off()

p1 <- ggplot(data_plot1, 
             aes(y = SNP, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(Bonferroni)')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/CHEERS/enrich_disease_individual_meta_bonferroni_immune.pdf', width = 8, height = 8)
p1
dev.off()

#####################
############
data_plot <- data_plot_bp[grep('IBD|UC|CD', data_plot_bp$SNP),]
data_plot <- data_plot[grep('All|Lead', data_plot$SNP),]
data_plot$snp <- 'All variants'
data_plot$snp[grep('Lead', data_plot$SNP)] <- 'Lead variants'
data_plot$snp <- factor(data_plot$snp, levels = c('Lead variants', 'All variants'))
data_plot$group <- 'IBD'
data_plot$group[grep('UC', data_plot$SNP)] <- 'UC'
data_plot$group[grep('CD', data_plot$SNP)] <- 'CD'
data_plot$snp_group <- paste0(data_plot$snp, ' (', data_plot$group, ')')

meta_sample = function(df, pheno){
    df0 = df[df$SNP == pheno,]

    meta_p = sapply(names(table(df0$celltype)), function(i){
        df00 = df0[df0$celltype == i,]
        meta_p = metap::sumlog(df00$p_value)
        return(meta_p$p)
    })
    
    table = data.frame(name = names(meta_p), p_value = meta_p, SNP = pheno, group = names(table(df0$group)), snp_group = names(table(df0$snp_group)))
    return(table)
}

meta_list = lapply(names(table(data_plot$SNP)), function(pheno){
    meta_sample(data_plot, pheno)
})

meta = Reduce(rbind, meta_list)

data_plot = meta

data_plot$log10_p_value <- -log10(data_plot$p_value)
hist(data_plot$log10_p_value)
data_plot$log10_p_value[data_plot$p_value >= 0.05] <- NA

data_plot1 <- data_plot[data_plot$group == 'IBD',]
data_plot1$snp_group <- factor(data_plot1$snp_group, levels = c('Lead variants (IBD)', 'All variants (IBD)'))
data_plot2 <- data_plot[data_plot$group == 'UC',]
data_plot2$snp_group <- factor(data_plot2$snp_group, levels = c('Lead variants (UC)', 'All variants (UC)'))
data_plot3 <- data_plot[data_plot$group == 'CD',]
data_plot3$snp_group <- factor(data_plot3$snp_group, levels = c('Lead variants (CD)', 'All variants (CD)'))

p1 <- ggplot(data_plot1, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)')+ 
  coord_fixed(expand = T)
p2 <- ggplot(data_plot2, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(P-value)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
p3 <- ggplot(data_plot3, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/CHEERS/enrich_disease_IBD_individual_meta.pdf', width = 10, height = 6.5)
p3/p2/p1
dev.off()

################
meta_list_fdr = meta_list
for (i in 1:6){
    meta_list_fdr[[i]]$p_val_adj = p.adjust(meta_list_fdr[[i]]$p_value, method = 'fdr')
}

meta = Reduce(rbind, meta_list_fdr)

data_plot = meta

data_plot$log10_p_value <- -log10(data_plot$p_val_adj)
hist(data_plot$log10_p_value)
data_plot$log10_p_value[data_plot$p_val_adj >= 0.05] <- NA

data_plot1 <- data_plot[data_plot$group == 'IBD',]
data_plot1$snp_group <- factor(data_plot1$snp_group, levels = c('Lead variants (IBD)', 'All variants (IBD)'))
data_plot2 <- data_plot[data_plot$group == 'UC',]
data_plot2$snp_group <- factor(data_plot2$snp_group, levels = c('Lead variants (UC)', 'All variants (UC)'))
data_plot3 <- data_plot[data_plot$group == 'CD',]
data_plot3$snp_group <- factor(data_plot3$snp_group, levels = c('Lead variants (CD)', 'All variants (CD)'))

p1 <- ggplot(data_plot1, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(FDR)')+ 
  coord_fixed(expand = T)
p2 <- ggplot(data_plot2, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(FDR)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
p3 <- ggplot(data_plot3, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(FDR)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/CHEERS/enrich_disease_IBD_individual_meta_fdr.pdf', width = 10, height = 6.5)
p3/p2/p1
dev.off()

################
meta_list_fdr = meta_list
for (i in 1:6){
    meta_list_fdr[[i]]$p_val_adj = p.adjust(meta_list_fdr[[i]]$p_value, method = 'bonferroni')
}

meta = Reduce(rbind, meta_list_fdr)
fwrite(meta, file = '../plots/SNP/CHEERS/enrich_disease_IBD_individual_meta_bonferroni.txt', col.names = T, row.names = F)

data_plot = meta

data_plot$log10_p_value <- -log10(data_plot$p_val_adj)
hist(data_plot$log10_p_value)
data_plot$log10_p_value[data_plot$p_val_adj >= 0.05] <- NA

data_plot1 <- data_plot[data_plot$group == 'IBD',]
data_plot1$snp_group <- factor(data_plot1$snp_group, levels = c('Lead variants (IBD)', 'All variants (IBD)'))
data_plot2 <- data_plot[data_plot$group == 'UC',]
data_plot2$snp_group <- factor(data_plot2$snp_group, levels = c('Lead variants (UC)', 'All variants (UC)'))
data_plot3 <- data_plot[data_plot$group == 'CD',]
data_plot3$snp_group <- factor(data_plot3$snp_group, levels = c('Lead variants (CD)', 'All variants (CD)'))

p1 <- ggplot(data_plot1, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(Bonferroni)')+ 
  coord_fixed(expand = T)
p2 <- ggplot(data_plot2, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(Bonferroni)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
p3 <- ggplot(data_plot3, 
             aes(y = snp_group, x = name, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,10), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(Bonferroni)')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/CHEERS/enrich_disease_IBD_individual_meta_bonferroni.pdf', width = 10, height = 6.5)
p3/p2/p1
dev.off()

pdf('../plots/SNP/CHEERS/enrich_disease_IBD_individual_meta_bonferroni_updated.pdf', width = 8, height = 6.5)
p3/(p2 + theme(legend.position = 'none'))/p1
dev.off()
