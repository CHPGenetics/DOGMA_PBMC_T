library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
library(stringr)
set.seed(1)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/')

addArchRThreads(threads = 64)

proj <- loadArchRProject('DOGMA_filtered_multiome/')

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

proj <- addMonocleTrajectory_fine('../monocle/cd4_naive.csv', 'MonocleCD4Naive_fine')
proj <- addMonocleTrajectory_fine('../monocle/cd4_memory.csv', 'MonocleCD4Memory_fine')
proj <- addMonocleTrajectory_fine('../monocle/cd8_naive.csv', 'MonocleCD8Naive_fine')
proj <- addMonocleTrajectory_fine('../monocle/cd8_memory.csv', 'MonocleCD8Memory_fine')
proj <- addMonocleTrajectory_fine('../monocle/mait.csv', 'MonocleMAIT_fine')
saveArchRProject(proj, outputDirectory = 'DOGMA_filtered_multiome')

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
  write.csv(regulator, file = paste0("DOGMA_filtered_multiome/Plots/Plot-", name, "-Traj-Heatmap_Cor.csv"), quote = F, col.names = F, row.names = F)
  
  ht1 <- plotTrajectoryHeatmap(trajGEM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
  ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder, labelTop = 100)
  plot <- ht1 + ht2
  plotPDF(plot, name = paste0("Plot-", name, "-Traj-Heatmap_Cor.pdf"), ArchRProj = proj, addDOC = FALSE, width = 12, height = 12)
}
monocle_analysis('MonocleCD4Naive_fine')
monocle_analysis('MonocleCD4Memory_fine')
monocle_analysis('MonocleCD8Naive_fine')
monocle_analysis('MonocleCD8Memory_fine')
monocle_analysis('MonocleMAIT_fine')

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
