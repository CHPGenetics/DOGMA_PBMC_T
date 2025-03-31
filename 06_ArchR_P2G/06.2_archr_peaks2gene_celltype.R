library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
set.seed(1)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/')

addArchRThreads(threads = 64)

proj <- loadArchRProject('DOGMA_filtered_multiome/')

p2g_analysis <- function(celltype){
  proj_sub <- proj[proj$celltype_id == celltype,]
  
  dir.create(file.path('DOGMA_filtered_multiome/', celltype))
  
  proj_sub <- addPeak2GeneLinks(
    ArchRProj = proj_sub,
    reducedDims = "Harmony_Combined",
    useMatrix = "GeneExpressionMatrix",
    maxDist = 1000000
  )
  
  dir.create(file.path('DOGMA_filtered_multiome/Peak2GeneLinks/', celltype))
  file.copy(from = 'DOGMA_filtered_multiome/Peak2GeneLinks/seRNA-Group-KNN.rds',
            to = paste0('DOGMA_filtered_multiome/Peak2GeneLinks/', celltype, '/seRNA-Group-KNN.rds'))
  file.copy(from = 'DOGMA_filtered_multiome/Peak2GeneLinks/seATAC-Group-KNN.rds',
            to = paste0('DOGMA_filtered_multiome/Peak2GeneLinks/', celltype, '/seATAC-Group-KNN.rds'))
  # saveArchRProject(proj, 'DOGMA_filtered_multiome')
  
  # p2g <- getPeak2GeneLinks(
  #   ArchRProj = proj_sub,
  #   corCutOff = 0.4,
  #   resolution = 1,
  #   returnLoops = FALSE
  # )
  # 
  # saveRDS(p2g, file = paste0('DOGMA_filtered_multiome/', celltype, '/p2g.RDS'))
  # 
  # p2g <- getPeak2GeneLinks(
  #   ArchRProj = proj_sub,
  #   corCutOff = 0.25,
  #   resolution = 1,
  #   returnLoops = FALSE
  # )
  # 
  # saveRDS(p2g, file = paste0('DOGMA_filtered_multiome/', celltype, '/p2g_cor_0.25.RDS'))
  # 
  # p2g <- getPeak2GeneLinks(
  #   ArchRProj = proj_sub,
  #   corCutOff = -1,
  #   FDRCutOff= 0.05,
  #   resolution = 1,
  #   returnLoops = FALSE
  # )
  # 
  # saveRDS(p2g, file = paste0('DOGMA_filtered_multiome/', celltype, '/p2g_fdr_0.05.RDS'))
  # 
  # p2g <- getPeak2GeneLinks(
  #   ArchRProj = proj_sub,
  #   corCutOff = -1,
  #   FDRCutOff= 0.05,
  #   resolution = 1
  # )
  # 
  # saveRDS(p2g, file = paste0('DOGMA_filtered_multiome/', celltype, '/p2g_fdr_0.05_loop.RDS'))
  # 
  # p2g <- getPeak2GeneLinks(
  #   ArchRProj = proj_sub,
  #   corCutOff = -1,
  #   FDRCutOff= 0.01,
  #   resolution = 1,
  #   returnLoops = FALSE
  # )
  # 
  # saveRDS(p2g, file = paste0('DOGMA_filtered_multiome/', celltype, '/p2g_fdr_0.01.RDS'))
  
  p2g <- getPeak2GeneLinks(
    ArchRProj = proj_sub,
    corCutOff = -1,
    FDRCutOff= 0.01,
    resolution = 1
  )
  
  saveRDS(p2g, file = paste0('DOGMA_filtered_multiome/', celltype, '/p2g_fdr_0.01_loop.RDS'))
  
  p2g <- getPeak2GeneLinks(
    ArchRProj = proj_sub,
    corCutOff = -1,
    FDRCutOff= 0.001,
    resolution = 1,
    returnLoops = FALSE
  )

  saveRDS(p2g, file = paste0('DOGMA_filtered_multiome/', celltype, '/p2g_fdr_0.001.RDS'))
  
  p2g <- getPeak2GeneLinks(
    ArchRProj = proj_sub,
    corCutOff = -1,
    FDRCutOff= 0.001,
    resolution = 1
  )
  
  saveRDS(p2g, file = paste0('DOGMA_filtered_multiome/', celltype, '/p2g_fdr_0.001_loop.RDS'))
}

for(i in names(table(proj$celltype_id))){
  p2g_analysis(i)
}

# celltype <- "06_Th17 cells"
# 
# proj_sub <- proj[proj$celltype_id == celltype,]
# 
# proj_sub <- addPeak2GeneLinks(
#   ArchRProj = proj_sub,
#   reducedDims = "Harmony_Combined",
#   useMatrix = "GeneExpressionMatrix",
#   maxDist = 1000000
# )
# 
# p2g <- getPeak2GeneLinks(
#   ArchRProj = proj_sub,
#   corCutOff = -1,
#   FDRCutOff= 0.05,
#   resolution = 1,
#   returnLoops = FALSE
# )
# 
# peak <- metadata(p2g)$peakSet
# gene <- metadata(p2g)$geneSet
# p2g_df <- as.data.frame(p2g)
# 
# peak_sub <- peak[p2g_df$idxATAC]
# gene_sub <- gene[p2g_df$idxRNA]
# peak_df <- data.frame(chr = seqnames(peak_sub), start = start(peak_sub), end = end(peak_sub))
# gene_df <- data.frame(chr = seqnames(gene_sub), pos = start(gene_sub), gene = gene_sub$name)
# peak_df$chr_pos <- paste(peak_df$chr, peak_df$start, peak_df$end, sep = '-')
# 
# p2g_df$peak_chr <- peak_df$chr
# p2g_df$peak_start <- peak_df$start
# p2g_df$peak_end <- peak_df$end
# p2g_df$peak_chr_pos <- peak_df$chr_pos
# p2g_df$gene_chr <- gene_df$chr
# p2g_df$gene_pos <- gene_df$pos
# p2g_df$gene <- gene_df$gene
# length(unique(p2g_df$peak_chr_pos))#216131
# length(unique(p2g_df$gene))#23946
# p2g_df$distance <- pmin(abs(p2g_df$peak_start - p2g_df$gene_pos), abs(p2g_df$peak_end - p2g_df$gene_pos))
# p2g_df$distance[p2g_df$peak_start <= p2g_df$gene_pos & p2g_df$gene_pos <= p2g_df$peak_end] <- 0
# 
# p2g_df <- p2g_df[order(p2g_df$FDR),]
# 
# rna <- readRDS("DOGMA_filtered_multiome/Peak2GeneLinks/seRNA-Group-KNN.rds")
# peak <- readRDS("DOGMA_filtered_multiome/Peak2GeneLinks/seATAC-Group-KNN.rds")
# 
# rna0 <- assays(rna)[[1]]
# peak0 <- assays(peak)[[1]]
# 
# plot_cor <- function(p2g){
#   rna0_0 <- rna0[p2g$idxRNA,]
#   peak0_0 <- peak0[p2g$idxATAC,]
#   rna_name <- p2g$gene
#   peak_name <- p2g$peak_chr_pos
#   
#   data_plot <- data.frame(rna = rna0_0, peak = peak0_0)
#   p <- ggplot(data_plot, aes(x = rna, y = peak))+
#     geom_point()+
#     geom_smooth(method = 'glm')+
#     stat_cor()+
#     theme_cowplot()+
#     labs(x = rna_name, y = peak_name)+
#     ggtitle('Peak-to-gene link in Th17 cells')+
#     theme(plot.title = element_text(hjust = 0.5))
#   return(p)
# }
# 
# plot_cor(p2g_df[which(p2g_df$Correlation < 0),][1,])
# p <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "celltype_id", corCutOff = 0,
#                           FDRCutOff= 0.05)
# plotPDF(p, name = "Peak2Gene-Heatmap", width = 16, height = 16, ArchRProj = proj, addDOC = FALSE)

