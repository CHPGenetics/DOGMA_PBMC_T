library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
set.seed(1)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/')

addArchRThreads(threads = 64)

proj <- loadArchRProject('DOGMA_filtered_multiome/')

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "celltype")

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered_multiome')

proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "celltype",
  pathToMacs2 = '/ihome/wchen/zhongli/.conda/envs/MACS2/bin/macs2'
)

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered_multiome')

getPeakSet(proj)

proj <- addPeakMatrix(proj)

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "JASPAR2020", name = "Motif", force = T)

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered_multiome')

proj <- addBgdPeaks(proj, force = T)

proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered_multiome')

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "celltype")

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGEM_MM <- correlateMatrices(proj, useMatrix1 = 'GeneExpressionMatrix', useMatrix2 = 'MotifMatrix', reducedDims = 'LSI_ATAC')

corGEM_MM$maxDelta <- rowData(seZ)[match(corGEM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGEM_MM <- corGEM_MM[order(abs(corGEM_MM$cor), decreasing = TRUE), ]
corGEM_MM <- corGEM_MM[which(!duplicated(gsub("\\-.*","",corGEM_MM[,"MotifMatrix_name"]))), ]
corGEM_MM$TFRegulator <- "No"
corGEM_MM$TFRegulator[which(corGEM_MM$cor > 0.5 & corGEM_MM$padj < 0.01 & corGEM_MM$maxDelta > quantile(corGEM_MM$maxDelta, 0.75))] <- "Positive"
corGEM_MM$TFRegulator[which(corGEM_MM$cor < -0.5 & corGEM_MM$padj < 0.01 & corGEM_MM$maxDelta > quantile(corGEM_MM$maxDelta, 0.75))] <- "Negative"
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


