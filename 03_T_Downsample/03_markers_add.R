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

plan('multicore')
options(future.globals.maxSize= 60*1024^3)

data <- readRDS('../output/tcell_annotated_updated.RDS')

Idents(data) <- 'predicted.celltype.l2'
p1 <- DimPlot(data, reduction = 'wnn2.umap', label = T, repel = T, label.size = 4) + guides(col = guide_legend(ncol = 1, override.aes = list(size=4)))
pdf('../plots/annotation_all_updated/predicted_legend.pdf', width = 10, height = 8)
print(p1)
dev.off()

data$Treg <- ifelse(data$celltype_updated %in% c('CD4+ Regulatory (Activated)', 'CD4+ Regulatory (Resting)'), 'CD4+ Regulatory', 'Other')
Idents(data) <- 'celltype_updated'
p1 <- FeatureScatter(object = data, feature1 = 'adt_CD4', feature2 = 'adt_CD25', pt.size = 0.2, plot.cor = F)+xlab('CD4 (ADT)')+ylab('CD25 (ADT)')
p2 <- FeatureScatter(object = data, feature1 = 'adt_CD4', feature2 = 'adt_CD25', pt.size = 0.2, plot.cor = F, group.by = 'Treg', cols = c('firebrick', 'lightgrey'))+xlab('CD4 (ADT)')+ylab('CD25 (ADT)')
pdf('../plots/annotation_all_updated/scatter_cd4_cd25.pdf', width = 18, height = 8)
print(p1|p2)
dev.off()

p1 <- FeatureScatter(object = data, feature1 = 'adt_CD4', feature2 = 'adt_CD127', pt.size = 0.2, plot.cor = F)+xlab('CD4 (ADT)')+ylab('CD127 (ADT)')
p2 <- FeatureScatter(object = data, feature1 = 'adt_CD4', feature2 = 'adt_CD127', pt.size = 0.2, plot.cor = F, group.by = 'Treg', cols = c('firebrick', 'lightgrey'))+xlab('CD4 (ADT)')+ylab('CD127 (ADT)')
pdf('../plots/annotation_all_updated/scatter_cd4_cd127.pdf', width = 18, height = 8)
print(p1|p2)
dev.off()

check_marker <- function(type, marker, name = str_split_fixed(marker, '-', 2)[,1]){
  if (type == 'rna'){
    p <- FeaturePlot(data, paste0('sct_', marker), reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle(paste0(name, ' (RNA)'))
  }
  if (type == 'adt'){
    p <- FeaturePlot(data, paste0('adt_', marker), reduction = 'wnn2.umap', cols = c("lightgrey", "darkgreen"), min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle(paste0(name, ' (ADT)'))
  }
  if (type == 'motif'){
    p <- FeaturePlot(data, paste0('chromvar_', marker), reduction = 'wnn2.umap', cols = c("lightgrey", "darkred"), min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle(paste0(name, ' (ATAC)'))
  }
  return(p)
}

pdf(paste0('../plots/annotation_all_updated/cd4_cd8.pdf'), width = 8, height = 4)
print(wrap_plots(check_marker('adt', 'CD4'), check_marker('adt', 'CD8'), ncol = 2))
dev.off()
