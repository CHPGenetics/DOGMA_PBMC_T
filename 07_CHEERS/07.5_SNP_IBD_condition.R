library(data.table)
library(stringr)
library(liftOver)
library(SNPlocs.Hsapiens.dbSNP141.GRCh38)
library(rtracklayer)
library(Signac)
library(cowplot)
library(ggsci)
library(ggpubr)
library(openxlsx)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')
# snp_all <- read.xlsx('../data/SNP/Table S1_Known IBD credible SNPs.xlsx', sheet = 1)
# snp_lead <- read.xlsx('../data/SNP/Table S1_Known IBD credible SNPs.xlsx', sheet = 2)
# snp_p <- read.xlsx('../data/SNP/Table S1_Known IBD credible SNPs.xlsx', sheet = 3)
# snp_p <- snp_p[-52,]
#
# huang <- read.xlsx('../data/SNP/41586_2017_BFnature22969_MOESM2_ESM.xlsx', sheet = 1)
# lange <- read.xlsx('../data/SNP/41588_2017_BFng3760_MOESM277_ESM.xlsx')
# lange_name <- as.character(lange[1,])
# lange <- lange[-c(1),]
# colnames(lange) <- lange_name
#
# huang <- huang[huang$tier2 =='No',]
# lange <- lange[lange$`Confidently fine-mapped in this study?` == T,]
#
# lange_ibd <- lange[lange$Phenotype == 'IBD',]$`CredibleSet (minimal SNP set with >95% posterior probability)`
# lange_ibd <- unlist(str_split(lange_ibd, ','))
# lange_uc <- lange[lange$Phenotype == 'UC',]$`CredibleSet (minimal SNP set with >95% posterior probability)`
# lange_uc <- unlist(str_split(lange_uc, ','))
# lange_cd <- lange[lange$Phenotype == 'CD',]$`CredibleSet (minimal SNP set with >95% posterior probability)`
# lange_cd <- unlist(str_split(lange_cd, ','))
# huang_ibd <- huang[huang$trait.reassigned == 'IBD',]$all.variant
# huang_ibd <- unlist(str_split(huang_ibd, ','))
# huang_uc <- huang[huang$trait.reassigned == 'UC',]$all.variant
# huang_uc <- unlist(str_split(huang_uc, ','))
# huang_cd <- huang[huang$trait.reassigned == 'CD',]$all.variant
# huang_cd <- unlist(str_split(huang_cd, ','))
#
# snp_ibd <- c(huang_ibd, lange_ibd)
# snp_uc <- c(huang_uc, lange_uc)
# snp_cd <- c(huang_cd, lange_cd)
#
# lift_snp <- function(snps){
#   snps <- snps[,c(2,3,3,1,4:ncol(snps))]
#   colnames(snps)[1:4] <- c('seqnames', 'start', 'end', 'SNP')
#   snps$trait_uc <- 0
#   snps$trait_cd <- 0
#   snps$trait_uc[snps$SNP %in% snp_uc] <- 1
#   snps$trait_cd[snps$SNP %in% snp_cd] <- 1
#
#   df <- snps[,c(1:4)]
#
#   cur <- makeGRangesFromDataFrame(df, keep.extra.columns = T)
#   ch = import.chain('../data/SNP/hg19ToHg38.over.chain')
#   ch
#
#   str(ch[[1]])
#
#   seqlevelsStyle(cur) = "UCSC"  # necessary
#   cur38 = liftOver(cur, ch)
#   class(cur38)
#
#   cur38 = unlist(cur38)
#   genome(cur38) = "hg38"
#   cur38
#
#   snps_38 <- data.frame(chr=seqnames(cur38),
#                         pos=start(cur38),
#                         SNP=cur38$SNP)
#
#   if(length(cur)-length(cur38)>0){
#     my_rsids <- setdiff(mcols(cur)$SNPS, mcols(cur38)$SNPS)
#
#     ## Note that the 1st call to snpsById() takes a long time but subsequent
#     ## calls are expected to be slightly faster.
#     my_snps <- snpsById(SNPlocs.Hsapiens.dbSNP141.GRCh38, my_rsids, ifnotfound="drop")
#     my_snps
#
#     snps_38 <- rbind(snps_38, data.frame(chr=gsub('ch', 'chr', seqnames(my_snps)),
#                                          pos=start(my_snps),
#                                          SNP=my_snps$RefSNP_id))
#   }
#
#
#   colnames(snps_38) <- c('chr_38', 'pos_38', 'SNP')
#   snps_38 <- snps_38[!duplicated(snps_38$SNP),]
#   snps <- merge(snps, snps_38, by = 'SNP', all.x = T)
#   snps$chr_38 <- as.character(snps$chr_38)
#   snps$chr_pos <- paste(snps$chr_38, snps$pos_38, sep = '-')
#   return(snps)
# }
#
# snp_all <- lift_snp(snp_all)
# snp_lead <- lift_snp(snp_lead)
# snp_p <- lift_snp(snp_p)
# 
# write.table(snp_all, '../output/SNP/SNP_IBD_all.txt', quote = F, row.names = F)
# write.table(snp_lead, '../output/SNP/SNP_IBD_lead.txt', quote = F, row.names = F)
# write.table(snp_p, '../output/SNP/SNP_IBD_p.txt', quote = F, row.names = F)

snp_all <- read.table('../output/SNP/SNP_IBD_all.txt', header = T)
snp_lead <- read.table('../output/SNP/SNP_IBD_lead.txt', header = T)
snp_p <- read.table('../output/SNP/SNP_IBD_p.txt', header = T)

snp_all <- snp_all[,-c(1:4)]
snp_all <- makeGRangesFromDataFrame(snp_all, seqnames.field = 'chr_38', start.field = 'pos_38', end.field = 'pos_38', keep.extra.columns = T)

snp_lead <- snp_lead[,-c(1:4)]
snp_lead <- makeGRangesFromDataFrame(snp_lead, seqnames.field = 'chr_38', start.field = 'pos_38', end.field = 'pos_38', keep.extra.columns = T)

snp_p <- snp_p[,-c(1:4)]
snp_p <- makeGRangesFromDataFrame(snp_p, seqnames.field = 'chr_38', start.field = 'pos_38', end.field = 'pos_38', keep.extra.columns = T)
###############
peak <- fread('../output/pseudo_bulk/DE_peak_DESeq2_paired.txt')

levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
           'CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
           'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1',
           'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17',
           'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh',
           'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other',
           'CD8+ Naive (Resting)', 'CD8+ Naive (Activated)',
           'CD8+ Regulatory',
           'CD8+ Memory (Resting)', 'CD8+ Memory (Activated)',
           'MAITs (Resting)', 'MAITs (Activated)', 'Gamma Delta'
)

df <- fread('../output/pseudo_bulk/peaks.annotation_sorted.txt')
df$name <- paste(df$Chr, df$Start-1, df$End, sep = '-')
i <- 1
anno <- df

result_combine_ps <- function(de, name){
  data <- de[de$cell_type == name,]
  print(table(data$gene == anno$name))
  data <- cbind(data, anno[,-c(1:4)])
  data <- data[order(data$p_val),]
  results = mutate(data, Expression = ifelse(data$p_val_adj < 0.05, "DEGs", "Not Significant"))
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$avg_logFC > 0) & (results$p_val_adj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$avg_logFC < 0) & (results$p_val_adj < 0.05))] <- "Down-regulated"
  # results$Expression[!results$gene %in% data2[data2$p_val_adj < 0.05,]$gene] <- "Not Significant"
  # p = ggplot(results, aes(avg_log2FC, -log10(p_val_adj))) + geom_point(aes(col = Expression)) +
  #   scale_color_manual(values = c("dodgerblue3", "darkgrey","firebrick")) +
  #   theme_bw(base_size = 12) +
  #   theme(text = element_text(size=16)) +
  #   theme(legend.position = "bottom") + geom_hline(yintercept = -log10(0.05), color = "grey1", lty = 2, lwd = 0.5)+
  #   ylab('-log10(FDR-adjusted P-value)')+
  #   xlab('log2(Fold Change)')+
  #   labs(col = '')+
  #   ggtitle(title)+
  #   theme(plot.title = element_text(hjust = 0.5))
  return(results)
}

peak_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  peak_list[[celltype]] <- result_combine_ps(peak, celltype)
}

markerList <- peak_list

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
  da_peaks <- markerList[[i]]
  da_peaks <- da_peaks[da_peaks$Expression != 'Not Significant',]
  marker_peak <- makeGRangesFromDataFrame(as.data.frame(str_split_fixed(da_peaks$gene, '-', 3)), 
                                          seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  # GCcontrol_peak <- makeGRangesFromDataFrame(as.data.frame(str_split_fixed(GCcontrolPeaks[[i]], '-', 3)), 
  #                                            seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  # Distcontrol_peak <- makeGRangesFromDataFrame(as.data.frame(str_split_fixed(DistcontrolPeaks[[i]], '-', 3)), 
  #                                              seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  data_plot <- rbind(data_plot, perform_fisher_test(marker_peak, snp_all, i))
  data_plot <- rbind(data_plot, perform_fisher_test(marker_peak, snp_lead, i))
  data_plot <- rbind(data_plot, perform_fisher_test(marker_peak, snp_p, i))
  # data_plot <- rbind(data_plot, perform_fisher_test(GCcontrol_peak, snp_all, i))
  # data_plot <- rbind(data_plot, perform_fisher_test(GCcontrol_peak, snp_lead, i))
  # data_plot <- rbind(data_plot, perform_fisher_test(GCcontrol_peak, snp_p, i))
  # data_plot <- rbind(data_plot, perform_fisher_test(Distcontrol_peak, snp_all, i))
  # data_plot <- rbind(data_plot, perform_fisher_test(Distcontrol_peak, snp_lead, i))
  # data_plot <- rbind(data_plot, perform_fisher_test(Distcontrol_peak, snp_p, i))
}
data_plot <- data_plot[-1,]
data_plot$snp <- gsub('snp_all', 'All variants', data_plot$snp)
data_plot$snp <- gsub('snp_lead', 'Leading variants', data_plot$snp)
data_plot$snp <- gsub('snp_p', 'Causal variants with > 95% certainty', data_plot$snp)
data_plot$peak <- gsub('marker_peak', 'Differentially accessible peaks', data_plot$peak)
# data_plot$peak <- gsub('GCcontrol_peak', 'GC-matched control peaks', data_plot$peak)
# data_plot$peak <- gsub('Distcontrol_peak', 'Distance-matched control peaks', data_plot$peak)
data_plot$snp <- factor(data_plot$snp, levels = c('Causal variants with > 95% certainty', 'Leading variants', 'All variants'))
# data_plot$peak <- factor(data_plot$peak, levels = c('Cell type-specific chromatin accessible peaks', 'GC-matched control peaks', 'Distance-matched control peaks'))

write.csv(data_plot, file = '../output/SNP_condition/enrich_immune_IBD.csv', col.names = T, row.names = F, quote = F)

# data_plot_immune <- data_plot[data_plot$snp == 'Causal variants with > 95% certainty',]
# data_plot_immune[which(data_plot_immune[data_plot_immune$peak == 'Cell type-specific chromatin accessible peaks',]$odds_ratio < data_plot_immune[data_plot_immune$peak == 'GC-matched control peaks',]$odds_ratio | 
#                          data_plot_immune[data_plot_immune$peak == 'Cell type-specific chromatin accessible peaks',]$odds_ratio < data_plot_immune[data_plot_immune$peak == 'Distance-matched control peaks',]$odds_ratio),]$name
#"02_CD4+ naive T cells (Resting)" "04_CD4+ naive T cells (Resting)" "05_Th1 cells" 
pdf('../plots/SNP_condition/enrich_immune_IBD.pdf', width = 8, height = 4)
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

pdf('../plots/SNP_condition/enrich_immune_paired_IBD.pdf', width = 10, height = 5)
ggplot(data_plot, aes(y = odds_ratio, x = peak, col = snp))+
  geom_violin(fill = 'white')+
  geom_point(fill = 'white', shape = 21, aes(size = -log10(p_value)))+
  geom_line(aes(group = name), col = 'darkgrey')+
  geom_hline(yintercept = 1, linetype = 'dashed', col = 'firebrick')+
  scale_color_npg()+
  coord_flip()+
  facet_wrap(~snp, nrow = 3, strip.position = 'left')+
  theme_cowplot()+
  theme(strip.text.y.left = element_blank(), strip.background = element_blank(), legend.position = 'bottom', legend.direction = 'vertical')+
  labs(x = 'Peaks', y = 'Odds ratio', col = 'SNPs', size = '-log10(P-value)')+
  scale_y_continuous(trans='log1p')
dev.off()

perform_fisher_test_split <- function(peak, snp, name, snp_name, group){
  overlap <- subsetByOverlaps(peak, snp)
  Overlap <- length(overlap)
  group2 <- length(snp)
  group1 <- length(peak)
  test <- fisher.test(matrix(c(Overlap, group2-Overlap, group1-Overlap, 15175044-group2-group1 +Overlap), 2, 2), alternative='greater')
  df <- data.frame(name = name, peak = deparse(substitute(peak)), snp = snp_name, odds_ratio = as.numeric(test$estimate), p_value = test$p.value, group = group)
  return(df)
}

data_plot <- as.data.frame(matrix(rep(NA, 7), ncol = 7))
colnames(data_plot) <- c('name', 'peak', 'snp', 'odds_ratio', 'p_value', 'group', 'order')
for(i in names(markerList)){
  da_peaks <- markerList[[i]]
  da_peaks <- da_peaks[da_peaks$Expression != 'Not Significant',]
  marker_peak <- makeGRangesFromDataFrame(as.data.frame(str_split_fixed(da_peaks$gene, '-', 3)), 
                                          seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  # GCcontrol_peak <- makeGRangesFromDataFrame(as.data.frame(str_split_fixed(GCcontrolPeaks[[i]], '-', 3)), 
  #                                            seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  # Distcontrol_peak <- makeGRangesFromDataFrame(as.data.frame(str_split_fixed(DistcontrolPeaks[[i]], '-', 3)), 
  #                                              seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  
  for(j in c('snp_all', 'snp_lead', 'snp_p')){
    disease_snps <- get(j)
    
    tmp <- perform_fisher_test_split(marker_peak, disease_snps, i, j, 'IBD')
    tmp$order <- order(tmp$odds_ratio, decreasing = T)
    data_plot <- rbind(data_plot, tmp)
    
    disease_snps <- get(j)
    disease_snps <- disease_snps[disease_snps$trait_uc == 1]
    
    tmp <- perform_fisher_test_split(marker_peak, disease_snps, i, j, 'UC')
    tmp$order <- order(tmp$odds_ratio, decreasing = T)
    data_plot <- rbind(data_plot, tmp)
    
    disease_snps <- get(j)
    disease_snps <- disease_snps[disease_snps$trait_cd == 1]
    
    tmp <- perform_fisher_test_split(marker_peak, disease_snps, i, j, 'CD')
    tmp$order <- order(tmp$odds_ratio, decreasing = T)
    data_plot <- rbind(data_plot, tmp)
  }
}
data_plot <- data_plot[-1,]
data_plot$snp <- gsub('snp_all', 'All variants', data_plot$snp)
data_plot$snp <- gsub('snp_lead', 'Leading variants', data_plot$snp)
data_plot$snp <- gsub('snp_p', 'Causal variants with > 95% certainty', data_plot$snp)
data_plot$peak <- gsub('marker_peak', 'Differentially accessible peaks', data_plot$peak)
# data_plot$peak <- gsub('GCcontrol_peak', 'GC-matched control peaks', data_plot$peak)
# data_plot$peak <- gsub('Distcontrol_peak', 'Distance-matched control peaks', data_plot$peak)
data_plot$snp <- factor(data_plot$snp, levels = c('Causal variants with > 95% certainty', 'Leading variants', 'All variants'))
# data_plot$peak <- factor(data_plot$peak, levels = c('Cell type-specific chromatin accessible peaks', 'GC-matched control peaks', 'Distance-matched control peaks'))
data_plot$snp_group <- paste0(data_plot$snp, ' (', data_plot$group, ')')

data_plot$log10_p_value <- -log10(data_plot$p_value)
data_plot$log10_p_value[data_plot$p_value >= 0.05] <- NA
data_plot$log10_p_value_control <- data_plot$log10_p_value
# data_plot$log10_p_value_control[which(data_plot$peak == 'Cell type-specific chromatin accessible peaks' & data_plot$order != 1)] <- NA
write.csv(data_plot, file = '../output/SNP_condition/enrich_disease_IBD.csv', col.names = T, row.names = F, quote = F)

data_plot$name <- factor(data_plot$name, levels = levels)
data_plot1 <- data_plot[data_plot$group == 'IBD',]
data_plot1$snp_group <- factor(data_plot1$snp_group, levels = c('Causal variants with > 95% certainty (IBD)', 'Leading variants (IBD)', 'All variants (IBD)'))
data_plot2 <- data_plot[data_plot$group == 'UC',]
data_plot2$snp_group <- factor(data_plot2$snp_group, levels = c('Causal variants with > 95% certainty (UC)', 'Leading variants (UC)', 'All variants (UC)'))
data_plot3 <- data_plot[data_plot$group == 'CD',]
data_plot3$snp_group <- factor(data_plot3$snp_group, levels = c('Causal variants with > 95% certainty (CD)', 'Leading variants (CD)', 'All variants (CD)'))

p1 <- ggplot(data_plot1[data_plot1$peak == 'Differentially accessible peaks',], 
             aes(y = snp_group, x = name, size = odds_ratio, fill = log10_p_value_control))+
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
p2 <- ggplot(data_plot2[data_plot2$peak == 'Differentially accessible peaks',], 
             aes(y = snp_group, x = name, size = odds_ratio, fill = log10_p_value_control))+
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
p3 <- ggplot(data_plot3[data_plot3$peak == 'Differentially accessible peaks',], 
             aes(y = snp_group, x = name, size = odds_ratio, fill = log10_p_value_control))+
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
pdf('../plots/SNP_condition/enrich_disease_IBD.pdf', width = 10, height = 7)
p3/p2/p1
dev.off()

#########
data_plot$snp <- gsub('Leading', 'Lead', data_plot$snp)
data_plot$snp_group <- gsub('Leading', 'Lead', data_plot$snp_group)
data_plot <- data_plot[data_plot$snp != 'Causal variants with > 95% certainty',]
data_plot1 <- data_plot[data_plot$group == 'IBD',]
data_plot1$snp_group <- factor(data_plot1$snp_group, levels = c('Lead variants (IBD)', 'All variants (IBD)'))
data_plot2 <- data_plot[data_plot$group == 'UC',]
data_plot2$snp_group <- factor(data_plot2$snp_group, levels = c('Lead variants (UC)', 'All variants (UC)'))
data_plot3 <- data_plot[data_plot$group == 'CD',]
data_plot3$snp_group <- factor(data_plot3$snp_group, levels = c('Lead variants (CD)', 'All variants (CD)'))

p1 <- ggplot(data_plot1[data_plot1$peak == 'Differentially accessible peaks',], 
             aes(y = snp_group, x = name, size = odds_ratio, fill = log10_p_value_control))+
  geom_point(shape = 21)+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,20), oob = scales::squish)+
  scale_size_area(
    max_size = 6,
    breaks = c(0,100,200),
    labels = c("0","100","200"),
    guide = "legend",
    limits = c(0, 200),
    oob = scales::squish
  )+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)', size = 'Odds ratio')+ 
  coord_fixed(expand = T)
p2 <- ggplot(data_plot2[data_plot2$peak == 'Differentially accessible peaks',], 
             aes(y = snp_group, x = name, size = odds_ratio, fill = log10_p_value_control))+
  geom_point(shape = 21)+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,20), oob = scales::squish)+
  scale_size_area(
    max_size = 6,
    breaks = c(0,100,200),
    labels = c("0","100","200"),
    guide = "legend",
    limits = c(0, 200),
    oob = scales::squish
  )+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(P-value)', size = 'Odds ratio')+ 
  rremove('x.text')+ 
  coord_fixed(expand = T)
p3 <- ggplot(data_plot3[data_plot3$peak == 'Differentially accessible peaks',], 
             aes(y = snp_group, x = name, size = odds_ratio, fill = log10_p_value_control))+
  geom_point(shape = 21)+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,20), oob = scales::squish)+
  scale_size_area(
    max_size = 6,
    breaks = c(0,100,200),
    labels = c("0","100","200"),
    guide = "legend",
    limits = c(0, 200),
    oob = scales::squish
  )+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')+
  labs(x = '', y = '', fill = '-log10(P-value)', size = 'Odds ratio')+ 
  rremove('x.text')+ 
  coord_fixed(expand = T)
pdf('../plots/SNP_condition/enrich_disease_IBD_updated.pdf', width = 10, height = 6.5)
p3/p2/p1
dev.off()

