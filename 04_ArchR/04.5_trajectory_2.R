library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
library(stringr)
set.seed(1)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/')

addArchRThreads(threads = 64)

proj <- loadArchRProject('DOGMA_filtered_multiome_2/')

###############################################
# Monocle3
###############################################
addMonocleTrajectory_fine <- function(path, name){
  table <- read.csv(path)
  monoclePT <- as.numeric(table$x)
  names(monoclePT) <- table$X
  names(monoclePT) <- gsub('_', '#', names(monoclePT))
  monoclePT <- ArchR:::.getQuantiles(monoclePT) * 100
  proj <- addCellColData(ArchRProj = proj, data = as.vector(monoclePT), 
                         name = name, cells = names(monoclePT), force = T)
  return(proj)
}

# #Add Monocle Myeloid Trajectory
# cds <- getMonocleTrajectories(
#   ArchRProj = proj,
#   groupBy = "celltype",
#   clusterParams = list(k = 50),
#   useGroups = c('Th1 cells',
#                 'Th17 cells', 'Tfh cells',
#                 'CD4+ memory T cells (Activated)', 'CD4+ memory T cells (Resting)',
#                 'Inducible regulatory T cells'),
#   principalGroup = "CD4+ memory T cells (Resting)",
#   embedding = "atac.umap"
# )
# 
# proj <- addMonocleTrajectory(
#   ArchRProj = proj,
#   monocleCDS = cds,
#   name = "MonocleCD4Memory",
#   groupBy = "celltype",
#   useGroups = c('Th1 cells',
#                 'Th17 cells', 'Tfh cells',
#                 'CD4+ memory T cells (Activated)', 'CD4+ memory T cells (Resting)',
#                 'Inducible regulatory T cells'),
#   force = TRUE
# )
# 
# trajMM  <- getTrajectory(ArchRProj = proj, name = "MonocleCD4Memory", useMatrix = "MotifMatrix", log2Norm = FALSE)
# 
# p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), labelRows = T, labelTop = 100)
# 
# plotPDF(p1, name = "Plot-MonocleCD4Memory-Traj-Heatmap.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)

proj <- addMonocleTrajectory_fine('../monocle_2/cd4_naive.csv', 'MonocleCD4Naive_fine')
proj <- addMonocleTrajectory_fine('../monocle_2/cd4_memory.csv', 'MonocleCD4Memory_fine')
proj <- addMonocleTrajectory_fine('../monocle_2/cd4_regulatory.csv', 'MonocleCD4Regulatory_fine')
proj <- addMonocleTrajectory_fine('../monocle_2/cd8_naive.csv', 'MonocleCD8Naive_fine')
proj <- addMonocleTrajectory_fine('../monocle_2/cd8_memory.csv', 'MonocleCD8Memory_fine')
proj <- addMonocleTrajectory_fine('../monocle_2/mait.csv', 'MonocleMAIT_fine')
saveArchRProject(proj, outputDirectory = 'DOGMA_filtered_multiome_2')

monocle_analysis <- function(name){
  trajMM  <- getTrajectory(ArchRProj = proj, name = name, useMatrix = "MotifMatrix", log2Norm = FALSE)
  p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), labelRows = T, labelTop = 75, varCutOff = 0.75)
  plotPDF(p1, name = paste0("Plot-", name, "-Traj-Heatmap.pdf"), ArchRProj = proj, addDOC = FALSE, width = 12, height = 12)
  trajGEM  <- getTrajectory(ArchRProj = proj, name = name, useMatrix = "GeneExpressionMatrix", log2Norm = FALSE)
  p1 <- plotTrajectoryHeatmap(trajGEM, pal = paletteContinuous(set = "solarExtra"), labelRows = T, labelTop = 75, varCutOff = 0.75)
  plotPDF(p1, name = paste0("Plot-", name, "-Traj-Heatmap_Expr.pdf"), ArchRProj = proj, addDOC = FALSE, width = 12, height = 12)
  trajPM  <- getTrajectory(ArchRProj = proj, name = name, useMatrix = "PeakMatrix", log2Norm = FALSE)
  p1 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"), labelRows = T, labelTop = 75, varCutOff = 0.75)
  plotPDF(p1, name = paste0("Plot-", name, "-Traj-Heatmap_Peak.pdf"), ArchRProj = proj, addDOC = FALSE, width = 12, height = 12)
  corGEM_MM <- correlateTrajectories(trajGEM, trajMM, corCutOff = 0.25, varCutOff1 = 0.65, varCutOff2 = 0.65)
  
  corGEM_MM[[1]]$matchname1
  corGEM_MM[[1]] <- corGEM_MM[[1]][str_split_fixed(corGEM_MM[[1]]$name2, ':', 2)[,1] == 'z',]
  
  trajGEM2 <- trajGEM[corGEM_MM[[1]]$name1, ]
  trajMM2 <- trajMM[corGEM_MM[[1]]$name2, ]
  
  trajCombined <- trajGEM2
  assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGEM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
  
  combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
  
  rowOrder <- match(rownames(combinedMat), rownames(trajGEM2))
  
  regulator <- rownames(trajGEM2)[rowOrder]
  regulator <- str_split_fixed(regulator, ':', 2)[,2]
  write.csv(regulator, file = paste0("DOGMA_filtered_multiome_2/Plots/Plot-", name, "-Traj-Heatmap_Cor.csv"), quote = F, col.names = F, row.names = F)
  
  ht1 <- plotTrajectoryHeatmap(trajGEM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
  ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
  plot <- ht1 + ht2
  plotPDF(plot, name = paste0("Plot-", name, "-Traj-Heatmap_Cor.pdf"), ArchRProj = proj, addDOC = FALSE, width = 12, height = 12)
}
monocle_analysis('MonocleCD4Naive_fine')
monocle_analysis('MonocleCD4Memory_fine')
monocle_analysis('MonocleCD4Regulatory_fine')
monocle_analysis('MonocleCD8Naive_fine')
monocle_analysis('MonocleCD8Memory_fine')
monocle_analysis('MonocleMAIT_fine')

proj <- addMonocleTrajectory_fine('../monocle_2/cd4_memory/Th17_Th1_1.csv', 'MonocleTh17_Th1_1')
proj <- addMonocleTrajectory_fine('../monocle_2/cd4_memory/Th17_Th1_2.csv', 'MonocleTh17_Th1_2')
proj <- addMonocleTrajectory_fine('../monocle_2/cd4_memory/Th17_Th1_3.csv', 'MonocleTh17_Th1_3')
proj <- addMonocleTrajectory_fine('../monocle_2/cd4_memory/Th17_Th1_4.csv', 'MonocleTh17_Th1_4')
proj <- addMonocleTrajectory_fine('../monocle_2/cd4_memory/Th17_Th1_5.csv', 'MonocleTh17_Th1_5')
proj <- addMonocleTrajectory_fine('../monocle_2/cd4_memory/Th17_Th1_6.csv', 'MonocleTh17_Th1_6')

monocle_analysis_cd4_memory <- function(name){
  trajMM  <- getTrajectory(ArchRProj = proj, name = name, useMatrix = "MotifMatrix", log2Norm = FALSE)
  p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), labelRows = T, labelTop = 75, varCutOff = 0.75)
  plotPDF(p1, name = paste0("cd4_memory/Plot-", name, "-Traj-Heatmap.pdf"), ArchRProj = proj, addDOC = FALSE, width = 12, height = 12)
  trajGEM  <- getTrajectory(ArchRProj = proj, name = name, useMatrix = "GeneExpressionMatrix", log2Norm = FALSE)
  p1 <- plotTrajectoryHeatmap(trajGEM, pal = paletteContinuous(set = "solarExtra"), labelRows = T, labelTop = 75, varCutOff = 0.75)
  plotPDF(p1, name = paste0("cd4_memory/Plot-", name, "-Traj-Heatmap_Expr.pdf"), ArchRProj = proj, addDOC = FALSE, width = 12, height = 12)
  trajPM  <- getTrajectory(ArchRProj = proj, name = name, useMatrix = "PeakMatrix", log2Norm = FALSE)
  p1 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"), labelRows = T, labelTop = 75, varCutOff = 0.75)
  plotPDF(p1, name = paste0("cd4_memory/Plot-", name, "-Traj-Heatmap_Peak.pdf"), ArchRProj = proj, addDOC = FALSE, width = 12, height = 12)
  corGEM_MM <- correlateTrajectories(trajGEM, trajMM, corCutOff = 0.25, varCutOff1 = 0.65, varCutOff2 = 0.65)
  
  corGEM_MM[[1]]$matchname1
  corGEM_MM[[1]] <- corGEM_MM[[1]][str_split_fixed(corGEM_MM[[1]]$name2, ':', 2)[,1] == 'z',]
  
  trajGEM2 <- trajGEM[corGEM_MM[[1]]$name1, ]
  trajMM2 <- trajMM[corGEM_MM[[1]]$name2, ]
  
  trajCombined <- trajGEM2
  assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGEM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
  
  combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
  
  rowOrder <- match(rownames(combinedMat), rownames(trajGEM2))
  
  regulator <- rownames(trajGEM2)[rowOrder]
  regulator <- str_split_fixed(regulator, ':', 2)[,2]
  write.csv(regulator, file = paste0("DOGMA_filtered_multiome_2/Plots/cd4_memory/Plot-", name, "-Traj-Heatmap_Cor.csv"), quote = F, col.names = F, row.names = F)
  
  ht1 <- plotTrajectoryHeatmap(trajGEM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
  ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
  plot <- ht1 + ht2
  plotPDF(plot, name = paste0("cd4_memory/Plot-", name, "-Traj-Heatmap_Cor.pdf"), ArchRProj = proj, addDOC = FALSE, width = 12, height = 12)
}
monocle_analysis_cd4_memory('MonocleTh17_Th1_1')
monocle_analysis_cd4_memory('MonocleTh17_Th1_2')
monocle_analysis_cd4_memory('MonocleTh17_Th1_3')
monocle_analysis_cd4_memory('MonocleTh17_Th1_4')
monocle_analysis_cd4_memory('MonocleTh17_Th1_5')
monocle_analysis_cd4_memory('MonocleTh17_Th1_6')

addTradeseqTrajectory_fine <- function(path, lineage, name){
  table <- read.csv(path, row.names = 'X')
  monoclePT <- as.numeric(table[,lineage])
  names(monoclePT) <- rownames(table)
  names(monoclePT) <- gsub('_', '#', names(monoclePT))
  monoclePT <- monoclePT[!is.na(monoclePT)]
  monoclePT <- ArchR:::.getQuantiles(monoclePT) * 100
  proj <- addCellColData(ArchRProj = proj, data = as.vector(monoclePT), 
                         name = name, cells = names(monoclePT), force = T)
  return(proj)
}
proj <- addTradeseqTrajectory_fine('../monocle_2/cd4_memory/tradeSeq_all.csv', 1, 'TradeSeq_1')
proj <- addTradeseqTrajectory_fine('../monocle_2/cd4_memory/tradeSeq_all.csv', 2, 'TradeSeq_2')
proj <- addTradeseqTrajectory_fine('../monocle_2/cd4_memory/tradeSeq_all.csv', 3, 'TradeSeq_3')
proj <- addTradeseqTrajectory_fine('../monocle_2/cd4_memory/tradeSeq_all.csv', 4, 'TradeSeq_4')

monocle_analysis_cd4_memory('TradeSeq_1')
monocle_analysis_cd4_memory('TradeSeq_2')
monocle_analysis_cd4_memory('TradeSeq_3')
monocle_analysis_cd4_memory('TradeSeq_4')

addTradeseqTrajectory_condition_fine <- function(path, lineage, condition, name){
  table <- read.csv(path, row.names = 'X')
  monoclePT <- as.numeric(table[,lineage])
  names(monoclePT) <- rownames(table)
  names(monoclePT) <- gsub('_', '#', names(monoclePT))
  monoclePT <- monoclePT[!is.na(monoclePT)]
  df <- as.data.frame(proj@cellColData)
  df <- df[names(monoclePT),]
  df_sub <- df[df$condition == condition,]
  monoclePT <- ArchR:::.getQuantiles(monoclePT) * 100
  monoclePT <- monoclePT[rownames(df_sub)]
  proj <- addCellColData(ArchRProj = proj, data = as.vector(monoclePT), 
                         name = name, cells = names(monoclePT), force = T)
  return(proj)
}
proj <- addTradeseqTrajectory_condition_fine('../monocle_2/cd4_memory/tradeSeq_all.csv', 1, 'Act_IL1B_IL23', 'TradeSeq_1_Act_IL1B_IL23')
proj <- addTradeseqTrajectory_condition_fine('../monocle_2/cd4_memory/tradeSeq_all.csv', 1, 'Act_IL1B_IL23_PGE2', 'TradeSeq_1_Act_IL1B_IL23_PGE2')

monocle_analysis_cd4_memory('TradeSeq_1_Act_IL1B_IL23')
monocle_analysis_cd4_memory('TradeSeq_1_Act_IL1B_IL23_PGE2')

addMonocleTrajectory_condition_fine <- function(path, condition, name){
  table <- read.csv(path)
  monoclePT <- as.numeric(table$x)
  names(monoclePT) <- table$X
  names(monoclePT) <- gsub('_', '#', names(monoclePT))
  df <- as.data.frame(proj@cellColData)
  df <- df[names(monoclePT),]
  df_sub <- df[df$condition == condition,]
  monoclePT <- ArchR:::.getQuantiles(monoclePT) * 100
  monoclePT <- monoclePT[rownames(df_sub)]
  proj <- addCellColData(ArchRProj = proj, data = as.vector(monoclePT), 
                         name = name, cells = names(monoclePT), force = T)
  return(proj)
}
proj <- addMonocleTrajectory_condition_fine('../monocle_2/cd4_memory/Th17_Th1_6.csv', 'Act_IL1B_IL23', 'MonocleTh17_Th1_6_Act_IL1B_IL23')
proj <- addMonocleTrajectory_condition_fine('../monocle_2/cd4_memory/Th17_Th1_6.csv', 'Act_IL1B_IL23_PGE2', 'MonocleTh17_Th1_6_Act_IL1B_IL23_PGE2')

monocle_analysis_cd4_memory_condition <- function(name, matrix){
  name1 <- paste0(name, '_Act_IL1B_IL23')
  name2 <- paste0(name, '_Act_IL1B_IL23_PGE2')
  
  trajMM1  <- getTrajectory(ArchRProj = proj, name = name1, useMatrix = matrix, log2Norm = FALSE)
  trajMM2  <- getTrajectory(ArchRProj = proj, name = name2, useMatrix = matrix, log2Norm = FALSE)
  
  p1 <- plotTrajectoryHeatmap(trajMM1, returnMat = TRUE, varCutOff = 0.85)
  p2 <- plotTrajectoryHeatmap(trajMM2, returnMat = TRUE, varCutOff = 0.85)
  table(rownames(p1) %in% rownames(p2))
  
  motif_share <- union(rownames(p1), rownames(p2))
  trajMM1_1 <- trajMM1[motif_share,]
  trajMM2_1 <- trajMM2[motif_share,]
  
  trajCombined <- trajMM1_1
  assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajMM1_1), 1, scale)) + t(apply(assay(trajMM2_1), 1, scale))
  
  combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
  
  rowOrder <- match(rownames(combinedMat), rownames(trajMM1_1))
  
  regulator <- rownames(trajMM1_1)[rowOrder]
  
  ht1 <- plotTrajectoryHeatmap(trajMM1_1,  pal = paletteContinuous(set = "solarExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
  ht1@ht_list[2] <- NULL
  ht2 <- plotTrajectoryHeatmap(trajMM2_1,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
  plot <- ht1 + ht2
  
  plotPDF(plot, name = paste0("cd4_memory/condition/Plot-", name, "-Traj-Heatmap_", matrix, ".pdf"), ArchRProj = proj, addDOC = FALSE, width = 14, height = 14)
}

monocle_analysis_cd4_memory_condition('MonocleTh17_Th1_6', 'MotifMatrix')
monocle_analysis_cd4_memory_condition('MonocleTh17_Th1_6', 'GeneExpressionMatrix')
monocle_analysis_cd4_memory_condition('MonocleTh17_Th1_6', 'PeakMatrix')

# proj <- addMonocleTrajectory(
#   ArchRProj = proj,
#   monocleCDS = cds,
#   name = "MonocleCD4Naive",
#   groupBy = "celltype",
#   useGroups = c("CD4+ naive T cells (Activated)", "CD4+ naive T cells (Resting)"),
#   force = TRUE
# )
# 
# head(proj$MonocleCD4Naive[!is.na(proj$MonocleCD4Naive)]) #NA means not in trajectory
# 
# trajMM  <- getTrajectory(ArchRProj = proj, name = "MonocleCD4Naive", useMatrix = "MotifMatrix", log2Norm = FALSE)
# 
# p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
# 
# plotPDF(p1, name = "Plot-MonocleCD4Naive-Traj-Heatmap.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)
# 
# proj <- addTrajectory(
#   ArchRProj = proj, 
#   name = "CD4Naive", 
#   groupBy = "celltype",
#   trajectory = c("CD4+ naive T cells (Activated)", "CD4+ naive T cells (Resting)"), 
#   embedding = "atac.umap", 
#   force = TRUE
# )
# 
# head(proj$CD4Naive[!is.na(proj$CD4Naive)]) #NA means not in trajectory
# 
# trajMM  <- getTrajectory(ArchRProj = proj, name = "CD4Naive", useMatrix = "MotifMatrix", log2Norm = FALSE)
# 
# p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
# 
# plotPDF(p1, name = "Plot-CD4Naive-Traj-Heatmap.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)
