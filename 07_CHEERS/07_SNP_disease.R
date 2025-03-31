library(data.table)
library(stringr)
library(liftOver)
library(SNPlocs.Hsapiens.dbSNP141.GRCh38)
library(rtracklayer)
library(Signac)
library(cowplot)
library(ggsci)
library(ggpubr)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')
snps <- readxl::read_xls('../data/SNP/41586_2015_BFnature13835_MOESM8_ESM.xls')

df <- snps[,c(4,5,5,3)]
colnames(df) <- c('seqnames', 'start', 'end', 'SNP')

cur <- makeGRangesFromDataFrame(df, keep.extra.columns = T)
ch = import.chain('../data/SNP/hg19ToHg38.over.chain')
ch

str(ch[[1]])

seqlevelsStyle(cur) = "UCSC"  # necessary
cur38 = liftOver(cur, ch)
class(cur38)

cur38 = unlist(cur38)
genome(cur38) = "hg38"
cur38

length(cur)-length(cur38)
my_rsids <- setdiff(mcols(cur)$SNP, mcols(cur38)$SNPs)

## Note that the 1st call to snpsById() takes a long time but subsequent
## calls are expected to be slightly faster.
my_snps <- snpsById(SNPlocs.Hsapiens.dbSNP141.GRCh38, my_rsids, ifnotfound="drop")
my_snps

snps_38 <- data.frame(chr=seqnames(cur38),
                      pos=start(cur38),
                      SNP=cur38$SNP)
snps_38 <- rbind(snps_38, data.frame(chr=gsub('ch', 'chr', seqnames(my_snps)),
                                     pos=start(my_snps),
                                     SNP=my_snps$RefSNP_id))
colnames(snps_38) <- c('chr_38', 'pos_38', 'SNP')
snps_38 <- snps_38[!duplicated(snps_38$SNP),]
snps <- merge(snps, snps_38, by = 'SNP', all.x = T)
snps$Group <- 'Immune'
snps$Group[snps$Disease %in% c('Platelet_counts', 'Red_blood_cell_traits',
                               'Urate_levels', 'Triglycerides', 'C_reactive_protein',
                               'HDL_cholesterol', 'LDL_cholesterol',
                               'Liver_enzyme_levels_gamma_glutamyl_transferase', 'Creatinine_levels')] <- 'Non-immune 1'
snps$Group[snps$Disease %in% c('Migraine', 'Renal_function_related_traits_BUN', 'Chronic_kidney_disease',
                               'Type_2_diabetes', 'Bone_mineral_density', 'Fasting_glucose_related_traits',
                               'Alzheimers_combined', 'Restless_legs_syndrome', 'Progressive_supranuclear_palsy')] <- 'Non-immune 2'
snps$Group_upper <- snps$Group
snps$Group_upper <- gsub(' 1', '', snps$Group_upper)
snps$Group_upper <- gsub(' 2', '', snps$Group_upper)
snps$chr_38 <- as.character(snps$chr_38)
snps$chr_pos <- paste(snps$chr_38, snps$pos_38, sep = '-')

write.csv(snps, '../output/SNP/SNP_disease.csv', quote = F)
###############
immune_snps <- snps[snps$Group_upper == 'Immune',]
immune_snps <- immune_snps[!duplicated(immune_snps$chr_pos),]
immune_snps <- makeGRangesFromDataFrame(immune_snps, seqnames.field = 'chr_38', start.field = 'pos_38', end.field = 'pos_38')

nonimmune_snps <- snps[snps$Group_upper == 'Non-immune',]
nonimmune_snps <- nonimmune_snps[!duplicated(nonimmune_snps$chr_pos),]
nonimmune_snps <- makeGRangesFromDataFrame(nonimmune_snps, seqnames.field = 'chr_38', start.field = 'pos_38', end.field = 'pos_38')

markerList <- readRDS('../output/ArchR/DOGMA_filtered_multiome/markerList.RDS')
GCcontrolPeaks <- readRDS('../output/ArchR/DOGMA_filtered_multiome/GCcontrolPeaks_chr.RDS')
DistcontrolPeaks <- readRDS('../output/ArchR/DOGMA_filtered_multiome/DistcontrolPeaks_chr.RDS')
markerList <- markerList[-9]

perform_fisher_test <- function(peak, snp, name){
  overlap <- subsetByOverlaps(peak, snp)
  Overlap <- length(overlap)
  group2 <- length(snp)
  group1 <- length(peak)
  test <- fisher.test(matrix(c(Overlap, group2-Overlap, group1-Overlap, 15175044-group2-group1 +Overlap), 2, 2), alternative='greater')
  df <- data.frame(name = name, peak = deparse(substitute(peak)), snp = deparse(substitute(snp)), odds_ratio = as.numeric(test$estimate), p_value = test$p.value)
  return(df)
}

data_plot <- as.data.frame(matrix(rep(NA, 5), ncol = 5))
colnames(data_plot) <- c('name', 'peak', 'snp', 'odds_ratio', 'p_value')
for(i in names(markerList)){
  marker_peak <- makeGRangesFromDataFrame(markerList[[i]])
  GCcontrol_peak <- makeGRangesFromDataFrame(as.data.frame(str_split_fixed(GCcontrolPeaks[[i]], '-', 3)), 
                                             seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  Distcontrol_peak <- makeGRangesFromDataFrame(as.data.frame(str_split_fixed(DistcontrolPeaks[[i]], '-', 3)), 
                                               seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  data_plot <- rbind(data_plot, perform_fisher_test(marker_peak, immune_snps, i))
  data_plot <- rbind(data_plot, perform_fisher_test(marker_peak, nonimmune_snps, i))
  data_plot <- rbind(data_plot, perform_fisher_test(GCcontrol_peak, immune_snps, i))
  data_plot <- rbind(data_plot, perform_fisher_test(GCcontrol_peak, nonimmune_snps, i))
  data_plot <- rbind(data_plot, perform_fisher_test(Distcontrol_peak, immune_snps, i))
  data_plot <- rbind(data_plot, perform_fisher_test(Distcontrol_peak, nonimmune_snps, i))
}
data_plot <- data_plot[-1,]
data_plot$snp <- gsub('nonimmune_snps', 'Non-immune', data_plot$snp)
data_plot$snp <- gsub('immune_snps', 'Immune', data_plot$snp)
data_plot$peak <- gsub('marker_peak', 'Cell type-specific chromatin accessible peaks', data_plot$peak)
data_plot$peak <- gsub('GCcontrol_peak', 'GC-matched control peaks', data_plot$peak)
data_plot$peak <- gsub('Distcontrol_peak', 'Distance-matched control peaks', data_plot$peak)
data_plot$snp <- factor(data_plot$snp, levels = c('Immune', 'Non-immune'))
data_plot$peak <- factor(data_plot$peak, levels = c('Cell type-specific chromatin accessible peaks', 'GC-matched control peaks', 'Distance-matched control peaks'))

write.csv(data_plot, file = '../output/SNP/enrich_immune.csv', col.names = T, row.names = F, quote = F)

data_plot_immune <- data_plot[data_plot$snp == 'Immune',]
data_plot_immune[which(data_plot_immune[data_plot_immune$peak == 'Cell type-specific chromatin accessible peaks',]$odds_ratio < data_plot_immune[data_plot_immune$peak == 'GC-matched control peaks',]$odds_ratio | 
                         data_plot_immune[data_plot_immune$peak == 'Cell type-specific chromatin accessible peaks',]$odds_ratio < data_plot_immune[data_plot_immune$peak == 'Distance-matched control peaks',]$odds_ratio),]$name
#[1] "06_Th17 cells"                      "07_Tfh cells"                       "08_CD4+ memory T cells (Activated)"
pdf('../plots/SNP/enrich_immune.pdf', width = 8, height = 4)
ggplot(data_plot, aes(y = odds_ratio, x = snp, fill = peak))+
  geom_violin()+
  geom_point(color = 'black', shape = 21, aes(size = -log10(p_value)))+
  geom_hline(yintercept = 1, linetype = 'dashed', col = 'darkgrey')+
  scale_fill_npg()+
  coord_flip()+
  facet_wrap(~peak, nrow = 3, strip.position = 'left')+
  theme_cowplot()+
  theme(strip.text.y.left = element_blank(), strip.background = element_blank(), legend.position = 'bottom', legend.direction = 'vertical')+
  labs(x = 'SNPs', y = 'Odds ratio', fill = 'Peaks', size = '-log10(P-value)')+
  scale_y_continuous(trans='log1p')
dev.off()

pdf('../plots/SNP/enrich_immune_paired.pdf', width = 10, height = 5)
ggplot(data_plot, aes(y = odds_ratio, x = peak, col = snp))+
  geom_violin(fill = 'white')+
  geom_point(fill = 'white', shape = 21, aes(size = -log10(p_value)))+
  geom_line(aes(group = name), col = 'darkgrey')+
  geom_hline(yintercept = 1, linetype = 'dashed', col = 'firebrick')+
  scale_color_manual(values = c('dodgerblue3', 'firebrick'))+
  coord_flip()+
  facet_wrap(~snp, nrow = 3, strip.position = 'left')+
  theme_cowplot()+
  theme(strip.text.y.left = element_blank(), strip.background = element_blank(), legend.position = 'bottom', legend.direction = 'vertical')+
  labs(x = 'Peaks', y = 'Odds ratio', col = 'SNPs', size = '-log10(P-value)')+
  scale_y_continuous(trans='log1p')
dev.off()

perform_fisher_test_split <- function(peak, snp, name, snp_name){
  overlap <- subsetByOverlaps(peak, snp)
  Overlap <- length(overlap)
  group2 <- length(snp)
  group1 <- length(peak)
  test <- fisher.test(matrix(c(Overlap, group2-Overlap, group1-Overlap, 15175044-group2-group1 +Overlap), 2, 2), alternative='greater')
  df <- data.frame(name = name, peak = deparse(substitute(peak)), snp = snp_name, odds_ratio = as.numeric(test$estimate), p_value = test$p.value)
  return(df)
}

data_plot <- as.data.frame(matrix(rep(NA, 6), ncol = 6))
colnames(data_plot) <- c('name', 'peak', 'snp', 'odds_ratio', 'p_value', 'order')
for(i in names(markerList)){
  marker_peak <- makeGRangesFromDataFrame(markerList[[i]])
  GCcontrol_peak <- makeGRangesFromDataFrame(as.data.frame(str_split_fixed(GCcontrolPeaks[[i]], '-', 3)), 
                                             seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  Distcontrol_peak <- makeGRangesFromDataFrame(as.data.frame(str_split_fixed(DistcontrolPeaks[[i]], '-', 3)), 
                                               seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  
  for(j in names(table(snps$Disease))){
    disease_snps <- snps[snps$Disease == j,]
    disease_snps <- disease_snps[!duplicated(disease_snps$chr_pos),]
    disease_snps <- makeGRangesFromDataFrame(disease_snps, seqnames.field = 'chr_38', start.field = 'pos_38', end.field = 'pos_38')
    
    tmp <- rbind(perform_fisher_test_split(marker_peak, disease_snps, i, j),
                 perform_fisher_test_split(GCcontrol_peak, disease_snps, i, j),
                 perform_fisher_test_split(Distcontrol_peak, disease_snps, i, j))
    tmp$order <- order(tmp$odds_ratio, decreasing = T)
    data_plot <- rbind(data_plot, tmp)
  }
}
data_plot <- data_plot[-1,]
data_plot$peak <- gsub('marker_peak', 'Cell type-specific chromatin accessible peaks', data_plot$peak)
data_plot$peak <- gsub('GCcontrol_peak', 'GC-matched control peaks', data_plot$peak)
data_plot$peak <- gsub('Distcontrol_peak', 'Distance-matched control peaks', data_plot$peak)
data_plot$Group <- 'Immune'
data_plot$Group[data_plot$snp %in% c('Platelet_counts', 'Red_blood_cell_traits',
                                     'Urate_levels', 'Triglycerides', 'C_reactive_protein',
                                     'HDL_cholesterol', 'LDL_cholesterol',
                                     'Liver_enzyme_levels_gamma_glutamyl_transferase', 'Creatinine_levels')] <- 'Non-immune 1'
data_plot$Group[data_plot$snp %in% c('Migraine', 'Renal_function_related_traits_BUN', 'Chronic_kidney_disease',
                                     'Type_2_diabetes', 'Bone_mineral_density', 'Fasting_glucose_related_traits',
                                     'Alzheimers_combined', 'Restless_legs_syndrome', 'Progressive_supranuclear_palsy')] <- 'Non-immune 2'
data_plot$log10_p_value <- -log10(data_plot$p_value)
data_plot$log10_p_value[data_plot$p_value >= 0.05] <- NA
data_plot$log10_p_value_control <- data_plot$log10_p_value
data_plot$log10_p_value_control[which(data_plot$peak == 'Cell type-specific chromatin accessible peaks' & data_plot$order != 1)] <- NA
write.csv(data_plot, file = '../output/SNP/enrich_disease.csv', col.names = T, row.names = F, quote = F)

data_plot1 <- data_plot[data_plot$Group == 'Immune',]
data_plot1$snp <- factor(data_plot1$snp, levels = names(table(data_plot1$snp)))
data_plot2 <- data_plot[data_plot$Group == 'Non-immune 1',]
data_plot2$snp <- factor(data_plot2$snp, levels = names(table(data_plot2$snp)))
data_plot3 <- data_plot[data_plot$Group == 'Non-immune 2',]
data_plot3$snp <- factor(data_plot3$snp, levels = names(table(data_plot3$snp)))

library(viridis)
p1 <- ggplot(data_plot1[data_plot1$peak == 'Cell type-specific chromatin accessible peaks',], 
             aes(y = snp, x = name, size = odds_ratio, fill = log10_p_value_control))+
  geom_point(shape = 21)+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,20), oob = scales::squish)+
  scale_size_area(
    max_size = 6,
    breaks = c(0,100,200,300),
    labels = c("0","100","200","300"),
    guide = "legend",
    limits = c(0, 300),
    oob = scales::squish
  )+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)', size = 'Odds ratio')+ 
  coord_fixed(expand = T)
p2 <- ggplot(data_plot2[data_plot2$peak == 'Cell type-specific chromatin accessible peaks',], 
             aes(y = snp, x = name, size = odds_ratio, fill = log10_p_value_control))+
  geom_point(shape = 21)+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,20), oob = scales::squish)+
  scale_size_area(
    max_size = 6,
    breaks = c(0,100,200,300),
    labels = c("0","100","200","300"),
    guide = "legend",
    limits = c(0, 300),
    oob = scales::squish
  )+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(P-value)', size = 'Odds ratio')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
p3 <- ggplot(data_plot3[data_plot3$peak == 'Cell type-specific chromatin accessible peaks',], 
             aes(y = snp, x = name, size = odds_ratio, fill = log10_p_value_control))+
  geom_point(shape = 21)+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,20), oob = scales::squish)+
  scale_size_area(
    max_size = 6,
    breaks = c(0,100,200,300),
    labels = c("0","100","200","300"),
    guide = "legend",
    limits = c(0, 300),
    oob = scales::squish
  )+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)', size = 'Odds ratio')+
  rremove('x.text')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/enrich_disease.pdf', width = 12, height = 12)
p3/p2/p1
dev.off()

p1 <- ggplot(data_plot1[data_plot1$peak == 'Cell type-specific chromatin accessible peaks',], 
             aes(y = snp, x = name, size = odds_ratio, fill = log10_p_value_control))+
  geom_point(shape = 21)+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,20), oob = scales::squish)+
  scale_size_area(
    max_size = 6,
    breaks = c(0,100,200,300),
    labels = c("0","100","200","300"),
    guide = "legend",
    limits = c(0, 300),
    oob = scales::squish
  )+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(P-value)', size = 'Odds ratio')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP/enrich_disease_immune.pdf', width = 8, height = 8)
p1
dev.off()
