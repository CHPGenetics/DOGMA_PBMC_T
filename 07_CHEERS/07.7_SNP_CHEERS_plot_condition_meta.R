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
library(patchwork)

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

table <- read.dir('../output/SNP_2/CHEERS/AllPeaks/Individual_condition_analysis/', '_disease_enrichment_pValues.txt')

data_plot_bp <- table[[1]]
for(i in 2:length(table)){
  data_plot_bp <- rbind(data_plot_bp, table[[i]])
}
data_plot_bp$celltype <- paste(str_split_fixed(data_plot_bp$name, '_', 3)[,1], str_split_fixed(data_plot_bp$name, '_', 3)[,2], sep = '_')
data_plot_bp$sample_condition <- str_split_fixed(data_plot_bp$name, '_', 3)[,3]
data_plot_bp$sample <- str_remove(data_plot_bp$sample_condition, 'Act_IL1B_IL23_PGE2_|Act_IL1B_IL23_')
data_plot_bp$sample <- as.factor(data_plot_bp$sample)
data_plot_bp$condition <- str_split_fixed(data_plot_bp$sample_condition, '_SB', 2)[,1]
data_plot_bp$condition <- factor(data_plot_bp$condition, levels = c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'))

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
data_plot$snp_group <- factor(as.character(data_plot$snp_group), levels = c('Lead variants (IBD)', 'All variants (IBD)',
                                                                            'Lead variants (UC)', 'All variants (UC)',
                                                                            'Lead variants (CD)', 'All variants (CD)'))

for(i in 1:6){
  name <- names(table(data_plot$snp_group))[i]
  p <- ggplot(data_plot[data_plot$snp_group == name,], aes(x = -log10(p_value), y = condition, fill = condition))+
           geom_boxplot()+
           geom_point(color = 'black', shape = 21)+
           geom_vline(xintercept = -log10(0.05), linetype = 'dashed', col = 'firebrick')+
    scale_fill_manual(values = c('firebrick', 'dodgerblue3'))
  p <- facet(p, facet.by = 'celltype', ncol = 1, strip.position = 'left')+
           theme_cowplot()+
           ggtitle(name)+
           theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5), strip.text.y.left = element_text(angle = 0), strip.background.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
           labs(y = '', x = '-log10(P-value)', fill = 'Conditions')
  assign(paste0('p',i), p)
}

pdf('../plots/SNP_2/CHEERS/enrich_disease_IBD_individual.pdf', width = 18, height = 18)
wrap_plots(p1,p2,p3,p4,p5,p6, nrow = 3)
dev.off()

compare_sample = function(df, pheno){
  df0 = df[df$SNP == pheno,]
  
  df0$celltype_condition = paste(df0$celltype, df0$condition, sep = ':')
  
  meta_p = sapply(seq(1,39,2), function(i){
    df01 = df0[df0$celltype_condition == names(table(df0$celltype_condition))[i],]
    df02 = df0[df0$celltype_condition == names(table(df0$celltype_condition))[i+1],]
    shared = intersect(df01$sample, df02$sample)
    df01 = df01[match(shared, df01$sample),]
    df02 = df02[match(shared, df02$sample),]

    table(df01$sample == df02$sample)

    test = t.test(qnorm(1-df01$p_value), qnorm(1-df02$p_value), paired = T)
    
    return(test$p.value)
  })
  
  table = data.frame(name = names(table(df0$celltype)), p_value = meta_p, SNP = pheno, group = names(table(df0$group)), snp_group = names(table(as.character(df0$snp_group))))
  return(table)
}

compare_list = lapply(names(table(data_plot$SNP)), function(pheno){
  compare_sample(data_plot, pheno)
})

compare = Reduce(rbind, compare_list)
fwrite(compare, file = '../plots/SNP_2/CHEERS/enrich_disease_IBD_condition_individual_compare.txt', col.names = T, row.names = F)

compare_sample = function(df, pheno){
  df0 = df[df$SNP == pheno,]
  
  df0$celltype_condition = paste(df0$celltype, df0$condition, sep = ':')
  
  meta_p = sapply(seq(1,39,2), function(i){
    df01 = df0[df0$celltype_condition == names(table(df0$celltype_condition))[i],]
    df02 = df0[df0$celltype_condition == names(table(df0$celltype_condition))[i+1],]
    # shared = intersect(df01$sample, df02$sample)
    # df01 = df01[match(shared, df01$sample),]
    # df02 = df02[match(shared, df02$sample),]

    # table(df01$sample == df02$sample)

    test = ks.test(df01$p_value, df02$p_value)
    
    return(test$p.value)
  })
  
  table = data.frame(name = names(table(df0$celltype)), p_value = meta_p, SNP = pheno, group = names(table(df0$group)), snp_group = names(table(as.character(df0$snp_group))))
  return(table)
}

compare_list = lapply(names(table(data_plot$SNP)), function(pheno){
  compare_sample(data_plot, pheno)
})

compare = Reduce(rbind, compare_list)
fwrite(compare, file = '../plots/SNP_2/CHEERS/enrich_disease_IBD_condition_individual_compare_ks.txt', col.names = T, row.names = F)

meta_sample = function(df, pheno){
  df0 = df[df$SNP == pheno,]
  
  df0$celltype_condition = paste(df0$celltype, df0$condition, sep = ':')
  
  meta_p = sapply(names(table(df0$celltype_condition)), function(i){
    df00 = df0[df0$celltype_condition == i,]
    meta_p = metap::sumlog(df00$p_value)
    return(meta_p$p)
  })
  
  table = data.frame(name = names(meta_p), p_value = meta_p, SNP = pheno, group = names(table(df0$group)), snp_group = names(table(as.character(df0$snp_group))))
  table$celltype = str_split_fixed(table$name, ':', 2)[,1]
  table$condition = str_split_fixed(table$name, ':', 2)[,2]
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
data_plot$snp_group <- factor(data_plot$snp_group, levels = c('Lead variants (CD)', 'All variants (CD)', 
'Lead variants (UC)', 'All variants (UC)',
'Lead variants (IBD)', 'All variants (IBD)'))
data_plot$condition = factor(data_plot$condition, levels = c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'))

p <- ggplot(data_plot, 
            aes(y = condition, x = celltype, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,20), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'right', plot.title = element_text(hjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(P-value)')+ 
  coord_fixed(expand = T)
p <- facet(p, facet.by = 'snp_group', ncol = 1, strip.position = 'left')+
  theme(strip.text.y.left = element_text(angle = 0), strip.background = element_blank())

pdf('../plots/SNP_2/CHEERS/enrich_disease_IBD_condition_individual_meta.pdf', width = 10, height = 7)
p
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
data_plot$snp_group <- factor(data_plot$snp_group, levels = c('Lead variants (CD)', 'All variants (CD)', 
'Lead variants (UC)', 'All variants (UC)',
'Lead variants (IBD)', 'All variants (IBD)'))
data_plot$condition = factor(data_plot$condition, levels = c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'))

p <- ggplot(data_plot, 
            aes(y = condition, x = celltype, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,20), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'right', plot.title = element_text(hjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(FDR)')+ 
  coord_fixed(expand = T)
p <- facet(p, facet.by = 'snp_group', ncol = 1, strip.position = 'left')+
  theme(strip.text.y.left = element_text(angle = 0), strip.background = element_blank())

pdf('../plots/SNP_2/CHEERS/enrich_disease_IBD_condition_individual_meta_fdr.pdf', width = 10, height = 7)
p
dev.off()

################
meta_list_fdr = meta_list
for (i in 1:6){
  meta_list_fdr[[i]]$p_val_adj = p.adjust(meta_list_fdr[[i]]$p_value, method = 'bonferroni')
}

meta = Reduce(rbind, meta_list_fdr)
fwrite(meta, file = '../plots/SNP_2/CHEERS/enrich_disease_IBD_condition_individual_meta_bonferroni.txt', col.names = T, row.names = F)

data_plot = meta

data_plot$log10_p_value <- -log10(data_plot$p_val_adj)
hist(data_plot$log10_p_value)
data_plot$log10_p_value[data_plot$p_val_adj >= 0.05] <- NA
data_plot$snp_group <- factor(data_plot$snp_group, levels = c('Lead variants (CD)', 'All variants (CD)', 
'Lead variants (UC)', 'All variants (UC)',
'Lead variants (IBD)', 'All variants (IBD)'))
data_plot$condition = factor(data_plot$condition, levels = c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'))

p <- ggplot(data_plot, 
            aes(y = condition, x = celltype, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,20), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'right', plot.title = element_text(hjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(Bonferroni)')+ 
  coord_fixed(expand = T)
p <- facet(p, facet.by = 'snp_group', ncol = 1, strip.position = 'left')+
  theme(strip.text.y.left = element_text(angle = 0), strip.background = element_blank())

pdf('../plots/SNP_2/CHEERS/enrich_disease_IBD_condition_individual_meta_bonferroni.pdf', width = 10, height = 7)
p
dev.off()

pdf('../plots/SNP_2/CHEERS/enrich_disease_IBD_condition_individual_meta_bonferroni_updated.pdf', width = 8, height = 6.5)
p + theme(legend.position = 'none')
dev.off()

compare = fread('../plots/SNP_2/CHEERS/enrich_disease_IBD_condition_individual_compare.txt')
data_plot$celltype_SNP = paste(data_plot$celltype, data_plot$SNP, sep = ':')
compare$celltype_SNP = paste(compare$name, compare$SNP, sep = ':')
compare_re = compare[rep(1:120, each = 2),]
table(data_plot$celltype_SNP == compare_re$celltype_SNP)
data_plot$compare_p_value = compare_re$p_value

p <- ggplot(data_plot, 
            aes(y = condition, x = celltype, fill = log10_p_value))+
  geom_tile(aes(color = ifelse(compare_p_value < 0.05, 'red', 'white'), width=0.85, height=0.85), size = 1)+
  scale_color_manual(values = c('red' = 'red', 'white' = NA), na.value = 'darkgrey')+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,20), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'right', plot.title = element_text(hjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(Bonferroni)')+ 
  coord_fixed(expand = T)
p <- facet(p, facet.by = 'snp_group', ncol = 1, strip.position = 'left')+
  theme(strip.text.y.left = element_text(angle = 0), strip.background = element_blank())

pdf('../plots/SNP_2/CHEERS/enrich_disease_IBD_condition_individual_meta_bonferroni_updated_compare.pdf', width = 8, height = 6.5)
p + theme(legend.position = 'none')
dev.off()

compare = fread('../plots/SNP_2/CHEERS/enrich_disease_IBD_condition_individual_compare_ks.txt')
data_plot$celltype_SNP = paste(data_plot$celltype, data_plot$SNP, sep = ':')
compare$celltype_SNP = paste(compare$name, compare$SNP, sep = ':')
compare_re = compare[rep(1:120, each = 2),]
table(data_plot$celltype_SNP == compare_re$celltype_SNP)
data_plot$compare_p_value = compare_re$p_value

p <- ggplot(data_plot, 
            aes(y = condition, x = celltype, fill = log10_p_value))+
  geom_tile(aes(color = ifelse(compare_p_value < 0.05, 'red', 'white'), width=0.85, height=0.85), size = 1)+
  scale_color_manual(values = c('red' = 'red', 'white' = NA), na.value = 'darkgrey')+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,20), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'right', plot.title = element_text(hjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(Bonferroni)')+ 
  coord_fixed(expand = T)
p <- facet(p, facet.by = 'snp_group', ncol = 1, strip.position = 'left')+
  theme(strip.text.y.left = element_text(angle = 0), strip.background = element_blank())

pdf('../plots/SNP_2/CHEERS/enrich_disease_IBD_condition_individual_meta_bonferroni_updated_compare_ks.pdf', width = 8, height = 6.5)
p + theme(legend.position = 'none')
dev.off()
