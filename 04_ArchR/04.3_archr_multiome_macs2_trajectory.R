library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
set.seed(1)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/')

addArchRThreads(threads = 64)

proj <- loadArchRProject('DOGMA_filtered_multiome/')

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "celltype")

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGEM_MM <- correlateMatrices(proj, useMatrix1 = 'GeneExpressionMatrix', useMatrix2 = 'MotifMatrix', reducedDims = 'Harmony_Combined', removeFromName1 = c("underscore"))

corGEM_MM$maxDelta <- rowData(seZ)[match(corGEM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGEM_MM <- corGEM_MM[order(abs(corGEM_MM$cor), decreasing = TRUE), ]
corGEM_MM <- corGEM_MM[which(!duplicated(gsub("\\-.*","",corGEM_MM[,"MotifMatrix_name"]))), ]
corGEM_MM$TFRegulator <- "No"
corGEM_MM$TFRegulator[which(corGEM_MM$cor > 0.25 & corGEM_MM$padj < 0.01 & corGEM_MM$maxDelta > quantile(corGEM_MM$maxDelta, 0.75))] <- "Positive"
corGEM_MM$TFRegulator[which(corGEM_MM$cor < -0.25 & corGEM_MM$padj < 0.01 & corGEM_MM$maxDelta > quantile(corGEM_MM$maxDelta, 0.75))] <- "Negative"
sort(corGEM_MM[corGEM_MM$TFRegulator=="Positive",1])
sort(corGEM_MM[corGEM_MM$TFRegulator=="Negative",1])
write.csv(corGEM_MM, file = 'DOGMA_filtered_multiome/corGEM_MM.csv')

p <- ggplot(data.frame(corGEM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() +
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") +
  scale_color_manual(values = c("No"="darkgrey", "Positive"="firebrick3", "Negative" = "dodgerblue3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, max(corGEM_MM$maxDelta)*1.05)
  )

pdf('DOGMA_filtered_multiome/Plots/Correlation-Motif-Deviation-TF-Expression.pdf', width = 4, height = 4)
p+theme(legend.position = 'none')
dev.off()

###################

data <- readRDS('../tcell_annotated_updated.RDS')
seuratUMAP <- function(name, seurat_name){
  proj@embeddings[[name]] <- proj@embeddings$UMAP_Combined
  table(rownames(proj@embeddings[[name]]$df) %in% gsub('_', '#', rownames(data@reductions[[seurat_name]]@cell.embeddings)))
  df <- data@reductions[[seurat_name]]@cell.embeddings
  rownames(df) <- gsub('_', '#', rownames(df))
  df <- df[rownames(proj@embeddings[[name]]$df),]
  colnames(df) <- colnames(proj@embeddings$UMAP_Combined$df)
  table(rownames(proj@embeddings[[name]]$df) == rownames(df))
  proj@embeddings[[name]]$df <- df
  return(proj)
}

proj <- seuratUMAP('wnn.umap', 'wnn.umap')
proj <- seuratUMAP('rna.umap', 'rna.umap')
proj <- seuratUMAP('adt.umap', 'adt.umap')
proj <- seuratUMAP('atac.umap', 'atac.umap')
proj <- seuratUMAP('wnn2.umap', 'wnn2.umap')
#Plot Embedding
p1 <- plotEmbedding(proj, name = "celltype", embedding = "wnn2.umap", size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj, name = "celltype", embedding = "wnn.umap", size = 1.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = "celltype", embedding = "rna.umap", size = 1.5, labelAsFactors=F, labelMeans=F)
p4 <- plotEmbedding(proj, name = "celltype", embedding = "adt.umap", size = 1.5, labelAsFactors=F, labelMeans=F)
p5 <- plotEmbedding(proj, name = "celltype", embedding = "atac.umap", size = 1.5, labelAsFactors=F, labelMeans=F)

#Save Plot
plotPDF(p1, p2, p3, p4, p5, name = "UMAP-Celltype", addDOC = FALSE)

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered_multiome')

#############
proj <- addImputeWeights(proj, reducedDims = 'Harmony_Combined')

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered_multiome')

library(ggplot2)
library(cowplot)
library(patchwork)
motif <- getMatrixFromProject(proj, useMatrix = 'MotifMatrix')
saveRDS(motif, file = 'DOGMA_filtered_multiome/motif.RDS')

rna <- getMatrixFromProject(proj, useMatrix = 'GeneExpressionMatrix')
saveRDS(rna, file = 'DOGMA_filtered_multiome/rna.RDS')

corGEM_MM <- read.csv('DOGMA_filtered_multiome/corGEM_MM.csv', row.names = 'X')

plot_magic <- function(marker, markerMotif){
  markerMotif_name <- str_split_fixed(markerMotif, '_', 2)[,1]
  markerMotif <- paste0("z:", markerMotif)
  p1 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneExpressionMatrix", 
    name = marker, 
    embedding = "wnn2.umap",
    quantCut = c(0.01, 0.99),
    imputeWeights = getImputeWeights(proj)
  ) + theme_cowplot() + ggtitle(paste0(marker, ' (RNA)')) + labs(x = 'wnn2UMAP_1', y = 'wnn2UMAP_2', fill = '') + theme(plot.title = element_text(hjust = 0.5), legend.position = 'right', legend.direction = 'vertical')
  p2 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "MotifMatrix", 
    name = markerMotif, 
    embedding = "wnn2.umap",
    quantCut = c(0.01, 0.99),
    imputeWeights = getImputeWeights(proj)
  ) + theme_cowplot() + ggtitle(paste0(markerMotif_name, ' (ATAC)')) + labs(x = 'wnn2UMAP_1', y = 'wnn2UMAP_2', fill = '') + theme(plot.title = element_text(hjust = 0.5), legend.position = 'right', legend.direction = 'vertical')
  plot <- p1 / p2
  pdf(paste0('DOGMA_filtered_multiome/Plots/plot_magic/', marker, '.pdf'), width = 4, height = 8)
  print(plot)
  dev.off()
}

plot_original <- function(marker, markerMotif){
  markerMotif_name <- str_split_fixed(markerMotif, '_', 2)[,1]
  markerMotif <- paste0("z:", markerMotif)
  p1 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneExpressionMatrix", 
    name = marker, 
    embedding = "wnn2.umap",
    quantCut = c(0.01, 0.99),
    imputeWeights = NULL
  ) + theme_cowplot() + ggtitle(paste0(marker, ' (RNA)')) + labs(x = 'wnn2UMAP_1', y = 'wnn2UMAP_2', fill = '') + theme(plot.title = element_text(hjust = 0.5), legend.position = 'right', legend.direction = 'vertical')
  p2 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "MotifMatrix", 
    name = markerMotif, 
    embedding = "wnn2.umap",
    quantCut = c(0.01, 0.99),
    imputeWeights = NULL
  ) + theme_cowplot() + ggtitle(paste0(markerMotif_name, ' (ATAC)')) + labs(x = 'wnn2UMAP_1', y = 'wnn2UMAP_2', fill = '') + theme(plot.title = element_text(hjust = 0.5), legend.position = 'right', legend.direction = 'vertical')
  plot <- p1 / p2
  pdf(paste0('DOGMA_filtered_multiome/Plots/plot_original/', marker, '.pdf'), width = 4, height = 8)
  print(plot)
  dev.off()
}

for(i in 1:nrow(corGEM_MM[corGEM_MM$TFRegulator != 'No',])){
  plot_magic(corGEM_MM[corGEM_MM$TFRegulator != 'No',1][i], corGEM_MM[corGEM_MM$TFRegulator != 'No',2][i])
  plot_original(corGEM_MM[corGEM_MM$TFRegulator != 'No',1][i], corGEM_MM[corGEM_MM$TFRegulator != 'No',2][i])
}

# plot_magic('RORA', 'RORA_9')
# plot_original('RORA', 'RORA_9')

plot_magic_split <- function(marker, markerMotif){
  markerMotif_name <- str_split_fixed(markerMotif, '_', 2)[,1]
  markerMotif <- paste0("z:", markerMotif)
  p1 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneExpressionMatrix", 
    name = marker, 
    embedding = "wnn2.umap",
    quantCut = c(0.01, 0.99),
    imputeWeights = getImputeWeights(proj)
  ) + theme_cowplot() + ggtitle(paste0(marker, ' (RNA)')) + labs(x = 'wnn2UMAP_1', y = 'wnn2UMAP_2', fill = '') + theme(plot.title = element_text(hjust = 0.5), legend.position = 'right', legend.direction = 'vertical')
  p2 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "MotifMatrix", 
    name = markerMotif, 
    embedding = "wnn2.umap",
    quantCut = c(0.01, 0.99),
    imputeWeights = getImputeWeights(proj)
  ) + theme_cowplot() + ggtitle(paste0(markerMotif_name, ' (ATAC)')) + labs(x = 'wnn2UMAP_1', y = 'wnn2UMAP_2', fill = '') + theme(plot.title = element_text(hjust = 0.5), legend.position = 'right', legend.direction = 'vertical')
  # plot <- p1 / p2
  # pdf(paste0('DOGMA_filtered_multiome_tcell/Plots/plot_magic/', marker, '.pdf'), width = 4, height = 8)
  # print(plot)
  # dev.off()
  plot <- list(p1, p2)
  return(plot)
}
plot1 <- plot_magic_split('NFKB1','NFKB1_198')
plot2 <- plot_magic_split('REL','REL_16')
plot3 <- plot_magic_split('JUND','JUND_586')

plot4 <- plot_magic_split('ETS1','ETS1_183')
plot5 <- plot_magic_split('ELF2','ELF2_376')
plot6 <- plot_magic_split('IKZF1','IKZF1_397')

plot7 <- plot_magic_split('LEF1','LEF1_188')
plot8 <- plot_magic_split('CEBPD','CEBPD_558')
plot9 <- plot_magic_split('RORC','RORC_351')

plot10 <- plot_magic_split('TCF7L2','TCF7L2_44')
plot11 <- plot_magic_split('ETS2','ETS2_377')
plot12 <- plot_magic_split('CREB1','CREB1_559')

pdf(paste0('DOGMA_filtered_multiome/Plots/plot_magic/', 'plots', '.pdf'), width = 36, height = 8)
wrap_plots(plot1[[1]],plot2[[1]],plot3[[1]],plot4[[1]],plot5[[1]],plot6[[1]],plot7[[1]],plot8[[1]],plot9[[1]],
           plot1[[2]],plot2[[2]],plot3[[2]],plot4[[2]],plot5[[2]],plot6[[2]],plot7[[2]],plot8[[2]],plot9[[2]], nrow = 2)
dev.off()

pdf(paste0('DOGMA_filtered_multiome/Plots/plot_magic/', 'plots_updated', '.pdf'), width = 48, height = 8)
wrap_plots(plot1[[1]],plot2[[1]],plot3[[1]],plot4[[1]],plot5[[1]],plot6[[1]],plot7[[1]],plot8[[1]],plot9[[1]],plot10[[1]],plot11[[1]],plot12[[1]],
           plot1[[2]],plot2[[2]],plot3[[2]],plot4[[2]],plot5[[2]],plot6[[2]],plot7[[2]],plot8[[2]],plot9[[2]],plot10[[2]],plot11[[2]],plot12[[2]], nrow = 2)
dev.off()

pdf(paste0('DOGMA_filtered_multiome/Plots/plot_magic/', 'plots_updated_updated', '.pdf'), width = 24, height = 16)
wrap_plots(plot1[[1]],plot2[[1]],plot3[[1]],plot4[[1]],plot5[[1]],plot6[[1]],
           plot1[[2]],plot2[[2]],plot3[[2]],plot4[[2]],plot5[[2]],plot6[[2]],
           plot7[[1]],plot8[[1]],plot9[[1]],plot10[[1]],plot11[[1]],plot12[[1]],
           plot7[[2]],plot8[[2]],plot9[[2]],plot10[[2]],plot11[[2]],plot12[[2]], nrow = 4)
dev.off()
