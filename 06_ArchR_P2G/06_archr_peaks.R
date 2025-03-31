library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
set.seed(1)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/')

addArchRThreads(threads = 64)

proj <- loadArchRProject('DOGMA_filtered_multiome/')

proj <- addPeakMatrix(proj, force = T)

proj$celltype_id <- proj$celltype
proj$celltype_id[proj$celltype_id == 'CD4+ Naive (Resting)'] <- '01_CD4+ Naive (Resting)'
proj$celltype_id[proj$celltype_id == 'CD4+ Naive (Activated)'] <- '02_CD4+ Naive (Activated)'
proj$celltype_id[proj$celltype_id == 'CD4+ Regulatory (Resting)'] <- '03_CD4+ Regulatory (Resting)'
proj$celltype_id[proj$celltype_id == 'CD4+ Regulatory (Activated)'] <- '04_CD4+ Regulatory (Activated)'
proj$celltype_id[proj$celltype_id == 'CD4+ Memory (Resting) - Th1'] <- '05_CD4+ Memory (Resting) - Th1'
proj$celltype_id[proj$celltype_id == 'CD4+ Memory (Activated) - Th1'] <- '06_CD4+ Memory (Activated) - Th1'
proj$celltype_id[proj$celltype_id == 'CD4+ Memory (Resting) - Th17'] <- '07_CD4+ Memory (Resting) - Th17'
proj$celltype_id[proj$celltype_id == 'CD4+ Memory (Activated) - Th17'] <- '08_CD4+ Memory (Activated) - Th17'
proj$celltype_id[proj$celltype_id == 'CD4+ Memory (Resting) - Tfh'] <- '09_CD4+ Memory (Resting) - Tfh'
proj$celltype_id[proj$celltype_id == 'CD4+ Memory (Activated) - Tfh'] <- '10_CD4+ Memory (Activated) - Tfh'
proj$celltype_id[proj$celltype_id == 'CD4+ Memory (Resting) - Other'] <- '11_CD4+ Memory (Resting) - Other'
proj$celltype_id[proj$celltype_id == 'CD4+ Memory (Activated) - Other'] <- '12_CD4+ Memory (Activated) - Other'
proj$celltype_id[proj$celltype_id == 'CD8+ Naive (Resting)'] <- '13_CD8+ Naive (Resting)'
proj$celltype_id[proj$celltype_id == 'CD8+ Naive (Activated)'] <- '14_CD8+ Naive (Activated)'
proj$celltype_id[proj$celltype_id == 'CD8+ Regulatory'] <- '15_CD8+ Regulatory'
proj$celltype_id[proj$celltype_id == 'CD8+ Memory (Resting)'] <- '16_CD8+ Memory (Resting)'
proj$celltype_id[proj$celltype_id == 'CD8+ Memory (Activated)'] <- '17_CD8+ Memory (Activated)'
proj$celltype_id[proj$celltype_id == 'MAITs (Resting)'] <- '18_MAITs (Resting)'
proj$celltype_id[proj$celltype_id == 'MAITs (Activated)'] <- '19_MAITs (Activated)'
proj$celltype_id[proj$celltype_id == 'Gamma Delta'] <- '20_Gamma Delta'

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "celltype_id",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveRDS(markersPeaks, file = 'DOGMA_filtered_multiome/markersPeaks.RDS')

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered_multiome')

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
markerList

saveRDS(markerList, file = 'DOGMA_filtered_multiome/markerList.RDS')

markersPeaks <- readRDS('DOGMA_filtered_multiome/markersPeaks.RDS')
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5",
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 16, height = 16, ArchRProj = proj, addDOC = FALSE)

markerList <- readRDS('DOGMA_filtered_multiome/markerList.RDS')

# celltype_proj <- proj[proj$celltype_id == name,]
celltype_proj <- proj
peaks_mat <- getMatrixFromProject(celltype_proj, useMatrix = 'PeakMatrix')
rowRanges(peaks_mat)$distToTSS <- celltype_proj@peakSet$distToTSS
rowRanges(peaks_mat)$GC <- celltype_proj@peakSet$GC
peaks_mat_df <- data.frame(chr=seqnames(peaks_mat),
                           start=start(peaks_mat),
                           end=end(peaks_mat),
                           distToTSS = rowRanges(peaks_mat)$distToTSS,
                           GC = rowRanges(peaks_mat)$GC)
rownames(peaks_mat_df) <- GRangesToString(rowRanges(peaks_mat), sep = c("-", "-"))

save(peaks_mat, file = 'DOGMA_filtered_multiome/PeakMatrix.RData')
# select_gc_control <- function(name){
#   peaks <- markerList[[name]]
#   peaks_gr <- peaks[,c(1,3,4,2)]
#   peaks_gr <- makeGRangesFromDataFrame(peaks_gr, keep.extra.columns = T)
#   
#   peaks_mat_sub <- subsetByOverlaps(peaks_mat, peaks_gr)
#   
#   peaks_mat_sub_df <- data.frame(chr=seqnames(peaks_mat_sub),
#                              start=start(peaks_mat_sub),
#                              end=end(peaks_mat_sub),
#                              distToTSS = rowRanges(peaks_mat_sub)$distToTSS,
#                              GC = rowRanges(peaks_mat_sub)$GC)
#   rownames(peaks_mat_sub_df) <- GRangesToString(rowRanges(peaks_mat_sub), sep = c("-", "-"))
#   
#   features.choose <- peaks_mat_df[-which(rownames(peaks_mat_df) %in% rownames(peaks_mat_sub_df)),]
#   
#   gc_control <- MatchRegionStats(
#     meta.feature = features.choose,
#     query.feature = peaks_mat_sub_df,
#     features.match = "GC",
#     n = nrow(peaks_mat_sub_df)
#   )
#   return(gc_control)
# }
# gc_control <- list()
# for(i in names(table(proj$celltype_id))){
#   gc_control[[i]] <- select_gc_control(i)
# }
# saveRDS(gc_control, file = 'DOGMA_filtered_multiome/GCcontrolPeaks.RDS')
# 
# select_dist_control <- function(name){
#   peaks <- markerList[[name]]
#   peaks_gr <- peaks[,c(1,3,4,2)]
#   peaks_gr <- makeGRangesFromDataFrame(peaks_gr, keep.extra.columns = T)
#   
#   peaks_mat_sub <- subsetByOverlaps(peaks_mat, peaks_gr)
#   
#   peaks_mat_sub_df <- data.frame(chr=seqnames(peaks_mat_sub),
#                                  start=start(peaks_mat_sub),
#                                  end=end(peaks_mat_sub),
#                                  distToTSS = rowRanges(peaks_mat_sub)$distToTSS,
#                                  GC = rowRanges(peaks_mat_sub)$GC)
#   rownames(peaks_mat_sub_df) <- GRangesToString(rowRanges(peaks_mat_sub), sep = c("-", "-"))
#   
#   features.choose <- peaks_mat_df[-which(rownames(peaks_mat_df) %in% rownames(peaks_mat_sub_df)),]
#   
#   dist_control <- MatchRegionStats(
#     meta.feature = features.choose,
#     query.feature = peaks_mat_sub_df,
#     features.match = "distToTSS",
#     n = nrow(peaks_mat_sub_df)
#   )
#   return(dist_control)
# }
# dist_control <- list()
# for(i in names(table(proj$celltype_id))){
#   dist_control[[i]] <- select_dist_control(i)
# }
# saveRDS(dist_control, file = 'DOGMA_filtered_multiome/DistcontrolPeaks.RDS')

#############
select_gc_control_chr <- function(name){
  peaks <- markerList[[name]]
  peaks_gr <- peaks[,c(1,3,4,2)]
  peaks_gr <- makeGRangesFromDataFrame(peaks_gr, keep.extra.columns = T)
  
  peaks_mat_sub <- subsetByOverlaps(peaks_mat, peaks_gr)
  
  peaks_mat_sub_df <- data.frame(chr=seqnames(peaks_mat_sub),
                                 start=start(peaks_mat_sub),
                                 end=end(peaks_mat_sub),
                                 distToTSS = rowRanges(peaks_mat_sub)$distToTSS,
                                 GC = rowRanges(peaks_mat_sub)$GC)
  rownames(peaks_mat_sub_df) <- GRangesToString(rowRanges(peaks_mat_sub), sep = c("-", "-"))
  
  peaks_mat_sub_df_chr <- peaks_mat_sub_df[peaks_mat_sub_df$chr == 'chr1',]
  features.choose <- peaks_mat_df[-which(rownames(peaks_mat_df) %in% rownames(peaks_mat_sub_df_chr)),]
  features.choose <- features.choose[features.choose$chr == 'chr1',]
  gc_control <- MatchRegionStats(
    meta.feature = features.choose,
    query.feature = peaks_mat_sub_df_chr,
    features.match = "GC",
    n = nrow(peaks_mat_sub_df_chr)
  )
  
  try(for(i in c(2:22, 'X')){
    chr_id <- paste0('chr',i)
    peaks_mat_sub_df_chr <- peaks_mat_sub_df[peaks_mat_sub_df$chr == chr_id,]
    features.choose <- peaks_mat_df[-which(rownames(peaks_mat_df) %in% rownames(peaks_mat_sub_df_chr)),]
    features.choose <- features.choose[features.choose$chr == chr_id,]
    gc_control <- c(gc_control, MatchRegionStats(
      meta.feature = features.choose,
      query.feature = peaks_mat_sub_df_chr,
      features.match = "GC",
      n = nrow(peaks_mat_sub_df_chr)
    ))
  })
  
  print(length(gc_control) == nrow(peaks_mat_sub_df))
  return(gc_control)
}
gc_control <- list()
for(i in names(table(proj$celltype_id))[-9]){
  gc_control[[i]] <- select_gc_control_chr(i)
}
saveRDS(gc_control, file = 'DOGMA_filtered_multiome/GCcontrolPeaks_chr.RDS')

select_dist_control_chr <- function(name){
  peaks <- markerList[[name]]
  peaks_gr <- peaks[,c(1,3,4,2)]
  peaks_gr <- makeGRangesFromDataFrame(peaks_gr, keep.extra.columns = T)
  
  peaks_mat_sub <- subsetByOverlaps(peaks_mat, peaks_gr)
  
  peaks_mat_sub_df <- data.frame(chr=seqnames(peaks_mat_sub),
                                 start=start(peaks_mat_sub),
                                 end=end(peaks_mat_sub),
                                 distToTSS = rowRanges(peaks_mat_sub)$distToTSS,
                                 GC = rowRanges(peaks_mat_sub)$GC)
  rownames(peaks_mat_sub_df) <- GRangesToString(rowRanges(peaks_mat_sub), sep = c("-", "-"))
  
  peaks_mat_sub_df_chr <- peaks_mat_sub_df[peaks_mat_sub_df$chr == 'chr1',]
  features.choose <- peaks_mat_df[-which(rownames(peaks_mat_df) %in% rownames(peaks_mat_sub_df_chr)),]
  features.choose <- features.choose[features.choose$chr == 'chr1',]
  dist_control <- MatchRegionStats(
    meta.feature = features.choose,
    query.feature = peaks_mat_sub_df,
    features.match = "distToTSS",
    n = nrow(peaks_mat_sub_df_chr)
  )
  
  try(for(i in c(2:22, 'X')){
    chr_id <- paste0('chr',i)
    peaks_mat_sub_df_chr <- peaks_mat_sub_df[peaks_mat_sub_df$chr == chr_id,]
    features.choose <- peaks_mat_df[-which(rownames(peaks_mat_df) %in% rownames(peaks_mat_sub_df_chr)),]
    features.choose <- features.choose[features.choose$chr == chr_id,]
    dist_control <- c(dist_control, MatchRegionStats(
      meta.feature = features.choose,
      query.feature = peaks_mat_sub_df,
      features.match = "distToTSS",
      n = nrow(peaks_mat_sub_df_chr)
    ))
  })
  
  print(length(dist_control) == nrow(peaks_mat_sub_df))
  return(dist_control)
}
dist_control <- list()
for(i in names(table(proj$celltype_id))[-9]){
  dist_control[[i]] <- select_dist_control_chr(i)
}
saveRDS(dist_control, file = 'DOGMA_filtered_multiome/DistcontrolPeaks_chr.RDS')