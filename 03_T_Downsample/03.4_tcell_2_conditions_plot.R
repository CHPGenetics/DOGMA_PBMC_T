library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(data.table)
library(Matrix)
library(dplyr)
library(stringr)
library(future)
library(harmony)
library(clustree)
library(openxlsx)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')

data <- read.csv('../output/tcell_annotated_updated.csv', row.names = 'X')
dim(data)
data <- data[data$condition %in% c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'),]
dim(data)

data$celltype <- factor(data$celltype_updated, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
                                                          'CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
                                                          'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1',
                                                          'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17',
                                                          'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh',
                                                          'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other',
                                                          'CD8+ Naive (Resting)', 'CD8+ Naive (Activated)',
                                                          'CD8+ Regulatory',
                                                          'CD8+ Memory (Resting)', 'CD8+ Memory (Activated)',
                                                          'MAITs (Resting)', 'MAITs (Activated)', 'Gamma Delta'
))

table <- read.xlsx('../../Single_Cell_Multimodal_Omics_Experiments_202200721.xlsx')
table <- table[table$Experiment_Name %in% data$orig.ident,]
data$sample <- factor(data$sample, levels = unique(table$Human_Subject_SB_Identifier))
data$condition <- factor(data$condition, levels = c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'))

library(ggplot2)
library(ggsci)
data_plot <- table(data$celltype, data$condition)
data_plot <- reshape2::melt(prop.table(data_plot, 2))
colnames(data_plot) <- c('Clusters', 'Conditions', 'Proportions')
data_plot$Clusters <- as.factor(data_plot$Clusters)
data_plot$Proportions <- data_plot$Proportions*100

p1 <- ggplot(data_plot, aes(x = Clusters, y = Proportions, fill = Conditions))+
  geom_col(position = 'dodge', col = 'white')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0), legend.position = 'bottom') + xlab('') + ylab('Proportion (%)') +
  scale_fill_npg() + labs(fill = '')
pdf('../plots/2_conditions/proportion_annotated_updated.pdf', width = 8, height = 4)
print(p1)
dev.off()

proportion <- function(data, name){
  data_plot <- table(data$celltype, data$condition)
  data_plot <- reshape2::melt(prop.table(data_plot, 2))
  colnames(data_plot) <- c('Clusters', 'Conditions', 'Proportions')
  data_plot$Clusters <- as.factor(data_plot$Clusters)
  data_plot$Proportions <- data_plot$Proportions*100
  data_plot$Donor <- name
  return(data_plot)
}

data_plot_list <- lapply(levels(data$sample), function(name){
  prop <- proportion(data[data$sample == name,], name)
  return(prop)
})
data_plot <- Reduce(rbind, data_plot_list)
data_plot <- rbind(data_plot[data_plot$Conditions == 'Act_IL1B_IL23',],
                   data_plot[data_plot$Conditions == 'Act_IL1B_IL23_PGE2',])

library(ggpubr)

myboxplot <- function(cluster){
  p <- ggboxplot(data_plot[data_plot$Clusters == cluster,], x = "Conditions", y = "Proportions",
                 fill = "Conditions", add = 'jitter')+ rremove('x.text')+
    stat_compare_means(paired = TRUE, method = 't.test', comparisons = list(c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2")))+
    theme(legend.position = 'none', axis.title.x = element_text(size = 16), plot.title = element_text(hjust = 0.5)) + xlab('') + ylab('Proportion (%)') + ggtitle(cluster) +
    scale_fill_npg() + labs(fill = '')
  return(p)
}

for(i in 1:20){
  assign(paste0('p',i), myboxplot(levels(data$celltype)[i]))
}
p21 <- myboxplot(levels(data$celltype)[20]) + theme(legend.position = 'bottom', legend.text = element_text(size=16), legend.title = element_text(size=16))

pdf('../plots/2_conditions/proportion_annotated_updated_pvalue_ttest.pdf', width = 20, height = 16)
wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20,
           ncol = 5)
dev.off()

pdf('../plots/2_conditions/proportion_annotated_updated_legend.pdf', width = 12, height = 5)
p21
dev.off()

myboxplot <- function(cluster){
  p <- ggboxplot(data_plot[data_plot$Clusters == cluster,], x = "Conditions", y = "Proportions",
                 fill = "Conditions", add = 'jitter')+ rremove('x.text')+
    stat_compare_means(paired = TRUE, method = 'wilcox.test', comparisons = list(c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2")))+
    theme(legend.position = 'none', axis.title.x = element_text(size = 16), plot.title = element_text(hjust = 0.5)) + xlab('') + ylab('Proportion (%)') + ggtitle(cluster) +
    scale_fill_npg() + labs(fill = '')
  return(p)
}

for(i in 1:20){
  assign(paste0('p',i), myboxplot(levels(data$celltype)[i]))
}

pdf('../plots/2_conditions/proportion_annotated_updated_pvalue_wilcox.pdf', width = 20, height = 16)
wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20,
           ncol = 5)
dev.off()

mypaired <- function(cluster){
  p <- ggpaired(rbind(data_plot[data_plot$Clusters == cluster,]), 
                x = "Conditions", y = "Proportions",
                fill = "Conditions", line.color = "darkgray", line.size = 0.4)+ rremove('x.text')+
    stat_compare_means(method = 't.test', paired = TRUE, comparisons = list(c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2")))+
    theme(legend.position = 'none', axis.title.x = element_text(size = 16), plot.title = element_text(hjust = 0.5)) + xlab('') + ylab('Proportion (%)') + ggtitle(cluster) +
    scale_fill_npg() + labs(fill = '')
  return(p)
}

data_plot$Conditions <- factor(data_plot$Conditions, levels = c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2"))
p <- ggpaired(data_plot, 
              x = "Conditions", y = "Proportions",
              fill = "Conditions", line.color = "darkgray", line.size = 0.4)+ rremove('x.text')
p <- facet(p, facet.by = 'Clusters', ncol = 23, strip.position = 'bottom')+
  stat_compare_means(method = 't.test', label = "p.signif", paired = TRUE, comparisons = list(c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2")))+
  theme(legend.position = 'bottom', strip.text.x = element_text(angle = 90), strip.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab('') + ylab('Proportion (%)') +
  scale_fill_npg() + labs(fill = '')

pdf('../plots/2_conditions/proportion_annotated_updated_pvalue_ttest_paired_merged.pdf', width = 18, height = 6)
p
dev.off()

for(i in 1:20){
  assign(paste0('p',i), mypaired(levels(data$celltype)[i]))
}
p21 <- mypaired(levels(data$celltype)[20]) + theme(legend.position = 'bottom', legend.text = element_text(size=16), legend.title = element_text(size=16))

pdf('../plots/2_conditions/proportion_annotated_updated_pvalue_ttest_paired.pdf', width = 20, height = 16)
wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20,
           ncol = 5)
dev.off()

pdf('../plots/2_conditions/proportion_annotated_updated_legend_paired.pdf', width = 12, height = 5)
p21
dev.off()

mypaired <- function(cluster){
  p <- ggpaired(rbind(data_plot[data_plot$Clusters == cluster,]), 
                x = "Conditions", y = "Proportions",
                fill = "Conditions", line.color = "darkgray", line.size = 0.4)+ rremove('x.text')+
    stat_compare_means(method = 'wilcox.test', paired = TRUE, comparisons = list(c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2")))+
    theme(legend.position = 'none', axis.title.x = element_text(size = 16), plot.title = element_text(hjust = 0.5)) + xlab('') + ylab('Proportion (%)') + ggtitle(cluster) +
    scale_fill_npg() + labs(fill = '')
  return(p)
}

p <- ggpaired(data_plot, 
              x = "Conditions", y = "Proportions",
              fill = "Conditions", line.color = "darkgray", line.size = 0.4)+ rremove('x.text')
p <- facet(p, facet.by = 'Clusters', ncol = 23, strip.position = 'bottom')+
  stat_compare_means(method = 'wilcox.test', label = "p.signif", paired = TRUE, comparisons = list(c("Act_IL1B_IL23", "Act_IL1B_IL23_PGE2")))+
  theme(legend.position = 'bottom', strip.text.x = element_text(angle = 90), strip.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab('') + ylab('Proportion (%)') +
  scale_fill_npg() + labs(fill = '')

pdf('../plots/2_conditions/proportion_annotated_updated_pvalue_wilcox_paired_merged.pdf', width = 18, height = 6)
p
dev.off()

for(i in 1:20){
  assign(paste0('p',i), mypaired(levels(data$celltype)[i]))
}

pdf('../plots/2_conditions/proportion_annotated_updated_pvalue_wilcox_paired.pdf', width = 20, height = 16)
wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20,
           ncol = 5)
dev.off()

###################
data <- readRDS('../output/tcell_annotated_updated_2_conditions.RDS')

library(openxlsx)
read.marker <- function(path, sheet){
  marker <- read.xlsx(xlsxFile = path, sheet = sheet)
  colnames(marker) <- marker[3,]
  marker <- marker[-c(1:3),]
  marker <- marker %>% mutate(across(contains('mean'), as.numeric)) %>% mutate(across(contains('LFC'), as.numeric)) %>% mutate(across(contains('p_'), as.numeric))
  return(marker)
}

th0.m.marker <- read.marker('../data/41467_2020_15543_MOESM6_ESM.xlsx', sheet = 3)
th0.m.marker.up <- th0.m.marker[th0.m.marker$LFC > 1,]
th0.m.marker.down <- th0.m.marker[th0.m.marker$LFC < -1,]

th0.n.marker <- read.marker('../data/41467_2020_15543_MOESM6_ESM.xlsx', sheet = 1)
th0.n.marker.up <- th0.n.marker[th0.n.marker$LFC > 1,]
th0.n.marker.down <- th0.n.marker[th0.n.marker$LFC < -1,]

th0.m.marker.5d <- read.marker('../data/41467_2020_15543_MOESM6_ESM.xlsx', sheet = 4)
th0.m.marker.5d.up <- th0.m.marker.5d[th0.m.marker.5d$LFC > 1,]
th0.m.marker.5d.down <- th0.m.marker.5d[th0.m.marker.5d$LFC < -1,]

th0.n.marker.5d <- read.marker('../data/41467_2020_15543_MOESM6_ESM.xlsx', sheet = 2)
th0.n.marker.5d.up <- th0.n.marker.5d[th0.n.marker.5d$LFC > 1,]
th0.n.marker.5d.down <- th0.n.marker.5d[th0.n.marker.5d$LFC < -1,]

th0.marker.up <- intersect(intersect(th0.m.marker.up$gene_name, th0.n.marker.up$gene_name), intersect(th0.m.marker.5d.up$gene_name, th0.n.marker.5d.up$gene_name))
th0.marker.down <- intersect(intersect(th0.m.marker.down$gene_name, th0.n.marker.down$gene_name), intersect(th0.m.marker.5d.down$gene_name, th0.n.marker.5d.down$gene_name))

DefaultAssay(data) <- "SCT"
data <- AddModuleScore(
  object = data,
  features = list(th0.marker.up),
  name = 'Activated'
)
data <- AddModuleScore(
  object = data,
  features = list(th0.marker.down),
  name = 'Resting'
)
p1 <- FeaturePlot(data, 'Activated1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + ggtitle('Activated')
p2 <- FeaturePlot(data, 'Resting1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99', raster = T) + ggtitle('Resting') 
pdf('../plots/2_conditions/activated_resting.pdf', width = 8, height = 4)
print(p1 | p2)
dev.off()

p3 <- DimPlot(data, reduction = 'wnn2.umap', group.by = 'condition', raster = T)
pdf('../plots/2_conditions/conditions.pdf', width = 6, height = 4)
p3
dev.off()

p4 <- DimPlot(data, reduction = 'wnn2.umap', group.by = 'condition', raster = T, cells = colnames(data)[data$condition %in% c("Act_IL1B_IL23")], cols = hue_pal()(2)[1]) + ggtitle('Act_IL1B_IL23') + theme(legend.position = 'none')
p5 <- DimPlot(data, reduction = 'wnn2.umap', group.by = 'condition', raster = T, cells = colnames(data)[data$condition %in% c("Act_IL1B_IL23_PGE2")], cols = hue_pal()(2)[2]) + ggtitle('Act_IL1B_IL23_PGE2') + theme(legend.position = 'none')
pdf('../plots/2_conditions/conditions_split.pdf', width = 8, height = 4)
print(p4 | p5)
dev.off()

pdf('../plots/2_conditions/activated_resting_condition.pdf', width = 11, height = 11)
wrap_plots(p1,p2,p4,p5, ncol = 2)
dev.off()

pdf('../plots/2_conditions/activated_resting_condition_updated.pdf', width = 15, height = 5)
wrap_plots(p1,p2,p3 + theme(legend.position = 'bottom'), ncol = 3)
dev.off()
