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

load('../output/tcell_downsample_predicted.RData')

plan('multisession')

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

read.marker <- function(path, sheet){
  marker <- read.xlsx(xlsxFile = path, sheet = sheet)
  colnames(marker) <- marker[3,]
  marker <- marker[-c(1:3),]
  marker <- marker %>% mutate(across(contains('mean'), as.numeric)) %>% mutate(across(contains('LFC'), as.numeric)) %>% mutate(across(contains('p_'), as.numeric))
  return(marker)
}

annotation <- function(data, dataset){
  
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
  
  DefaultAssay(data) <- 'SCT'
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
  p1 <- FeaturePlot(data, 'Activated1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle('Activated')
  p2 <- FeaturePlot(data, 'Resting1', reduction = 'wnn2.umap', min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle('Resting') 
  pdf(paste0('../plots/', dataset, '/activated_resting.pdf'), width = 8, height = 4)
  print(p1 | p2)
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/naive_memory.pdf'), width = 8, height = 8)
  print(wrap_plots(check_marker('adt', 'CD4'), check_marker('adt', 'CD8'),
                   check_marker('adt', 'CD45RA'), check_marker('adt', 'CD45RO'), ncol = 2))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/CD4_naive.pdf'), width = 8, height = 8)
  print(wrap_plots(check_marker('adt', 'CD4'), check_marker('adt', 'CD45RA'), 
                   check_marker('adt', 'CD62L'), check_marker('rna', 'CCR7'), ncol = 2))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/CD8_naive.pdf'), width = 8, height = 8)
  print(wrap_plots(check_marker('adt', 'CD8'), check_marker('adt', 'CD45RA'), 
                   check_marker('adt', 'CD62L'), check_marker('rna', 'CCR7'), ncol = 2))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/activated_Th.pdf'), width = 16, height = 8)
  print(wrap_plots(check_marker('adt', 'CD4'), check_marker('adt', 'CD44'), 
                   check_marker('adt', 'CD25'), check_marker('adt', 'CD69'),
                   check_marker('adt', 'CD71'), 
                   check_marker('rna', 'SLC3A2'), check_marker('rna', 'SLC7A5'), ncol = 4))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/Treg.pdf'), width = 16, height = 16)
  print(wrap_plots(check_marker('adt', 'CD152'), 
                   check_marker('adt', 'CD279'), check_marker('adt', 'CD39'),
                   check_marker('adt', 'CD103'),
                   check_marker('adt', 'CD25'), check_marker('adt', 'CD134'),
                   check_marker('adt', 'CD137'), check_marker('adt', 'CD223'), 
                   check_marker('rna', 'CCR7'), check_marker('motif', 'MA0850.1', 'FOXP3'), 
                   check_marker('rna', 'FOXP3'), check_marker('rna', 'IKZF2'), 
                   check_marker('rna', 'SMAD3'), check_marker('rna', 'AHR'), 
                   check_marker('rna', 'CTLA4'), check_marker('rna', 'IL2RA'), ncol = 4))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/Th17.pdf'), width = 20, height = 16)
  print(wrap_plots(check_marker('adt', 'CD161'), check_marker('adt', 'CD278'),
                   check_marker('adt', 'CD194'), check_marker('adt', 'CD196'),
                   check_marker('rna', 'IL21R'), check_marker('rna', 'IL23R'), 
                   check_marker('rna', 'IL17F'), 
                   check_marker('rna', 'TNF'), check_marker('rna', 'CCL20'), 
                   check_marker('motif', 'MA1151.1', 'RORC'),
                   check_marker('motif', 'MA0071.1', 'RORA'), check_marker('motif', 'MA0072.1', 'RORA'), check_marker('rna', 'RORA'), 
                   check_marker('motif', 'MA1634.1', 'BATF'), check_marker('rna', 'BATF'), 
                   check_marker('motif', 'MA1419.1', 'IRF4'), check_marker('rna', 'IRF4'), 
                   check_marker('motif', 'MA1419.1', 'IRF4'), check_marker('rna', 'IRF4'), ncol = 5))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/Th1.pdf'), width = 16, height = 12)
  print(wrap_plots(check_marker('adt', 'CD183'), check_marker('adt', 'CD195'),
                   check_marker('adt', 'CD26'), check_marker('adt', 'CD94'),
                   check_marker('adt', 'CD278'), 
                   check_marker('rna', 'IL18R1'), check_marker('rna', 'IFNG'), 
                   check_marker('rna', 'LTA'), check_marker('rna', 'TNF'), 
                   check_marker('motif', 'MA0690.1', 'TBX21'), check_marker('rna', 'TBX21'), ncol = 4))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/Tfh.pdf'), width = 20, height = 16)
  print(wrap_plots(check_marker('adt', 'CD185'), check_marker('adt', 'CD84'),
                   check_marker('adt', 'CD150'), check_marker('adt', 'CD200'),
                   check_marker('adt', 'CD272'), check_marker('adt', 'CD278'), 
                   check_marker('adt', 'CD279'), check_marker('adt', 'TIGIT'), 
                   check_marker('adt', 'CD57'), check_marker('adt', 'CD154'), 
                   check_marker('adt', 'CD304'),
                   check_marker('rna', 'IL6R'), check_marker('rna', 'TNFSF8'), 
                   check_marker('rna', 'IL21R'), check_marker('rna', 'CD40LG'), 
                   check_marker('motif', 'MA0463.2', 'BCL6'), check_marker('rna', 'BCL6'), ncol = 5))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/Mono.pdf'), width = 16, height = 12)
  print(wrap_plots(check_marker('adt', 'CD14'), check_marker('adt', 'CD16'),
                   check_marker('adt', 'CD192'), check_marker('adt', 'CD83'),
                   check_marker('adt', 'CD123'), check_marker('adt', 'CD99'), 
                   check_marker('adt', 'CD64'), check_marker('adt', 'CD35'), 
                   check_marker('adt', 'CD11b'), check_marker('adt', 'CD88'), 
                   check_marker('adt', 'CX3CR1'),
                   check_marker('rna', 'CD83'), ncol = 4))
  dev.off()
}

annotation(data, 'annotation')

tcell <- data
rm(data)
DefaultAssay(tcell) <- 'SCT'

analysis <- function(data, dataset, cluster){
  dir.create(file.path(paste0('../output/', dataset, '/tcell/'), cluster), showWarnings = FALSE)
  dir.create(file.path(paste0('../plots/', dataset, '/tcell/'), cluster), showWarnings = FALSE)
  
  DefaultAssay(data) <- "RNA"
  data <- SCTransform(data, method = "glmGamPoi", vars.to.regress = "percent.mt") %>% RunPCA(reduction.name = "pca") %>%
    RunHarmony(group.by.vars = c("sample"), max.iter.harmony = 40, reduction = 'pca', assay.use = 'SCT',project.dim = FALSE,  reduction.save = "harmony_RNA")
  data <- RunUMAP(data, reduction = "harmony_RNA", dims = 1:30, assay = 'SCT', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
  data <- FindNeighbors(data, dims = 1:30, reduction = 'harmony_RNA')
  data <- FindClusters(data, resolution = 0.2, graph.name = 'SCT_snn', algorithm = 3)
  data$rna_clusters <- Idents(data)
  
  DefaultAssay(data) <- 'ADT'
  # we will use all ADT features for dimensional reduction
  # we set a dimensional reduction name to avoid overwriting the
  VariableFeatures(data) <- rownames(data[["ADT"]])
  data <- NormalizeData(data, normalization.method = 'CLR', margin = 2) %>%
    ScaleData() %>% RunPCA(reduction.name = 'apca') %>%
    RunHarmony(group.by.vars = c("sample", "cryopreservation"), max.iter.harmony = 40, reduction = 'apca', assay.use = 'ADT',project.dim = FALSE,  reduction.save = "harmony_ADT")
  data <- RunUMAP(data, reduction = "harmony_ADT", dims = 1:30, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
  data <- FindNeighbors(data, dims = 1:30, reduction = 'harmony_ADT')
  data <- FindClusters(data, resolution = 0.2, graph.name = 'ADT_snn', algorithm = 3)
  data$adt_clusters <- Idents(data)
  
  DefaultAssay(data) <- "peaks"
  data <- RunTFIDF(data)
  data <- FindTopFeatures(data, min.cutoff = 'q0')
  data <- RunSVD(data)
  data <- RunHarmony(data, group.by.vars = c("sample"), max.iter.harmony = 40, reduction = 'lsi', assay.use = 'peaks',project.dim = FALSE,  reduction.save = "harmony_peaks")
  data <- RunUMAP(data, reduction = "harmony_peaks", dims = 2:30, assay = 'peaks', reduction.name = 'atac.umap', reduction.key = 'atacUMAP_')
  data <- FindNeighbors(data, dims = 2:30, reduction = 'harmony_peaks')
  data <- FindClusters(data, resolution = 0.2, graph.name = 'peaks_snn', algorithm = 3)
  data$atac_clusters <- Idents(data)
  
  data <- FindMultiModalNeighbors(object = data,
                                  reduction.list = list("harmony_RNA", "harmony_peaks", "harmony_ADT"),
                                  dims.list = list(1:30, 2:30, 1:30))
  data <- RunUMAP(data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_" )
  data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = 0.2)
  data$wnn_clusters <- Idents(data)
  
  data <- FindMultiModalNeighbors(object = data,
                                  reduction.list = list("harmony_peaks", "harmony_ADT"),
                                  dims.list = list(2:30, 1:30),
                                  knn.graph.name = "wknn2",
                                  snn.graph.name = "wsnn2",
                                  weighted.nn.name = "weighted.nn2")
  data <- RunUMAP(data, nn.name = "weighted.nn2", reduction.name = "wnn2.umap", reduction.key = "wnn2UMAP_" )
  data <- FindClusters(data, graph.name = "wsnn2", algorithm = 3, resolution = 0.2)
  data$wnn2_clusters <- Idents(data)
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/UMAP_condition_split.pdf'), width = 32, height = 40)
  p1 <- DimPlot(data, reduction = 'wnn.umap', split.by = 'condition', label = F, repel = TRUE, label.size = 4, ncol = 4)
  p2 <- DimPlot(data, reduction = 'rna.umap', split.by = 'condition', label = F, repel = TRUE, label.size = 4, ncol = 4)
  p3 <- DimPlot(data, reduction = 'adt.umap', split.by = 'condition', label = F, repel = TRUE, label.size = 4, ncol = 4)
  p4 <- DimPlot(data, reduction = 'atac.umap', split.by = 'condition', label = F, repel = TRUE, label.size = 4, ncol = 4)
  p5 <- DimPlot(data, reduction = 'wnn2.umap', split.by = 'condition', label = F, repel = TRUE, label.size = 4, ncol = 4)
  print(wrap_plots(p1, p2, p3, p4, p5, ncol = 1))
  dev.off()
  
  data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = c(seq(0.02, 0.2, 0.02)))
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/tree_WNN.pdf'), width = 10, height = 10)
  print(clustree(data@meta.data, prefix = 'wsnn_res.'))
  dev.off()
  
  data <- FindClusters(data, graph.name = "SCT_snn", algorithm = 3, resolution = c(seq(0.02, 0.2, 0.02)))
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/tree_RNA.pdf'), width = 10, height = 10)
  print(clustree(data@meta.data, prefix = 'SCT_snn_res.'))
  dev.off()
  
  data <- FindClusters(data, graph.name = "ADT_snn", algorithm = 3, resolution = c(seq(0.02, 0.2, 0.02)))
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/tree_ADT.pdf'), width = 10, height = 10)
  print(clustree(data@meta.data, prefix = 'ADT_snn_res.'))
  dev.off()
  
  data <- FindClusters(data, graph.name = "peaks_snn", algorithm = 3, resolution = c(seq(0.02, 0.2, 0.02)))
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/tree_ATAC.pdf'), width = 10, height = 10)
  print(clustree(data@meta.data, prefix = 'peaks_snn_res.'))
  dev.off()
  
  data <- FindClusters(data, graph.name = "wsnn2", algorithm = 3, resolution = c(seq(0.02, 0.2, 0.02)))
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/tree_WNN2.pdf'), width = 10, height = 10)
  print(clustree(data@meta.data, prefix = 'wsnn2_res.'))
  dev.off()
  
  p1 <- FeaturePlot(data, 'percent.mt', reduction = 'wnn.umap')
  p2 <- FeaturePlot(data, 'percent.mt', reduction = 'rna.umap')
  p3 <- FeaturePlot(data, 'percent.mt', reduction = 'adt.umap')
  p4 <- FeaturePlot(data, 'percent.mt', reduction = 'atac.umap')
  p5 <- FeaturePlot(data, 'percent.mt', reduction = 'wnn2.umap')
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/mito.pdf'), width = 24, height = 16)
  print(wrap_plots(p1, p2, p3, p4, p5, ncol = 3))
  dev.off()
  
  p1 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = 'wnn_clusters')
  p2 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/annotated_predicted_WNN.pdf'), width = 20, height = 8)
  print(p1 | p2)
  dev.off()
  
  p1 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = 'rna_clusters')
  p2 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/annotated_predicted_RNA.pdf'), width = 20, height = 8)
  print(p1 | p2)
  dev.off()
  
  p1 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = 'adt_clusters')
  p2 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/annotated_predicted_ADT.pdf'), width = 20, height = 8)
  print(p1 | p2)
  dev.off()
  
  p1 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = 'atac_clusters')
  p2 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/annotated_predicted_ATAC.pdf'), width = 20, height = 8)
  print(p1 | p2)
  dev.off()
  
  p1 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = 'wnn2_clusters')
  p2 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/annotated_predicted_WNN2.pdf'), width = 20, height = 8)
  print(p1 | p2)
  dev.off()
  
  saveRDS(data, file = paste0('../output/', dataset, '/tcell_', cluster, '.RDS'))
}

check_marker_cluster <- function(reduction, type, marker, name = str_split_fixed(marker, '-', 2)[,1]){
  if (type == 'rna'){
    p <- FeaturePlot(data, paste0('sct_', marker), reduction = reduction, min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle(paste0(name, ' (RNA)'))
  }
  if (type == 'adt'){
    p <- FeaturePlot(data, paste0('adt_', marker), reduction = reduction, cols = c("lightgrey", "darkgreen"), min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle(paste0(name, ' (ADT)'))
  }
  if (type == 'motif'){
    p <- FeaturePlot(data, paste0('chromvar_', marker), reduction = reduction, cols = c("lightgrey", "darkred"), min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle(paste0(name, ' (ATAC)'))
  }
  return(p)
}

annotation_cluster <- function(graph, threshold, reduction, dataset, cluster){
  library(openxlsx)
  Idents(data) <- paste0(graph, threshold)
  data@active.ident <- factor(data@active.ident, levels = 0:(length(levels(Idents(data)))-1))
  
  rna.markers <- FindAllMarkers(data, assay = 'SCT', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.xlsx(rna.markers, paste0('../output/', dataset, '/tcell/', cluster, '/rna_markers.xlsx'), colNames = T, rowNames = F, overwrite = T)
  rna.top10 <- rna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/heatmap_RNA.pdf'), width = 16, height = 16)
  print(DoHeatmap(data, assay = 'SCT', features = rna.top10$gene) + NoLegend())
  dev.off()
  
  adt.markers <- FindAllMarkers(data, assay = 'ADT', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.xlsx(adt.markers, paste0('../output/', dataset, '/tcell/', cluster, '/adt_markers.xlsx'), colNames = T, rowNames = F, overwrite = T)
  adt.top10 <- adt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/heatmap_ADT.pdf'), width = 16, height = 16)
  print(DoHeatmap(data, assay = 'ADT', features = adt.top10$gene) + NoLegend())
  dev.off()
  
  motif.markers <- FindAllMarkers(data, assay = 'chromvar', only.pos = TRUE,
                                  mean.fxn = rowMeans,
                                  fc.name = "avg_diff")
  motif.markers$motif <- ConvertMotifID(data, id = motif.markers$gene, assay = 'peaks')
  motif.markers$name <- paste0(motif.markers$motif, ' (', motif.markers$gene, ')')
  write.xlsx(motif.markers, paste0('../output/', dataset, '/tcell/', cluster, '/motif_markers.xlsx'), colNames = T, rowNames = F, overwrite = T)
  motif.top10 <- motif.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_diff)
  
  motif.markers_uni <- motif.markers[!duplicated(motif.markers$name),]
  motif <- data@assays$chromvar@data[motif.markers_uni$gene,]
  table(rownames(motif) == motif.markers_uni$gene)
  rownames(motif) <- motif.markers_uni$name
  data[['chromvar2']] <- CreateAssayObject(data = motif)
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/heatmap_motif.pdf'), width = 16, height = 16)
  print(DoHeatmap(data, assay = 'chromvar2', slot = 'data', features = motif.top10$name) + NoLegend())
  dev.off()
  
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
  
  DefaultAssay(data) <- 'SCT'
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
  p1 <- FeaturePlot(data, 'Activated1', reduction = reduction, min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle('Activated')
  p2 <- FeaturePlot(data, 'Resting1', reduction = reduction, min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle('Resting') 
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/activated_resting.pdf'), width = 8, height = 4)
  print(p1 | p2)
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/naive_memory.pdf'), width = 8, height = 8)
  try(print(wrap_plots(check_marker_cluster(reduction, 'adt', 'CD4'), check_marker_cluster(reduction, 'adt', 'CD8'),
                       check_marker_cluster(reduction, 'adt', 'CD45RA'), check_marker_cluster(reduction, 'adt', 'CD45RO'), ncol = 2)))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/CD4_naive.pdf'), width = 8, height = 8)
  try(print(wrap_plots(check_marker_cluster(reduction, 'adt', 'CD4'), check_marker_cluster(reduction, 'adt', 'CD45RA'), 
                       check_marker_cluster(reduction, 'adt', 'CD62L'), check_marker_cluster(reduction, 'rna', 'CCR7'), ncol = 2)))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/CD8_naive.pdf'), width = 8, height = 8)
  try(print(wrap_plots(check_marker_cluster(reduction, 'adt', 'CD8'), check_marker_cluster(reduction, 'adt', 'CD45RA'), 
                       check_marker_cluster(reduction, 'adt', 'CD62L'), check_marker_cluster(reduction, 'rna', 'CCR7'), ncol = 2)))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/activated_Th.pdf'), width = 16, height = 8)
  try(print(wrap_plots(check_marker_cluster(reduction, 'adt', 'CD4'), check_marker_cluster(reduction, 'adt', 'CD44'), 
                       check_marker_cluster(reduction, 'adt', 'CD25'), check_marker_cluster(reduction, 'adt', 'CD69'),
                       check_marker_cluster(reduction, 'adt', 'CD71'), 
                       check_marker_cluster(reduction, 'rna', 'SLC3A2'), check_marker_cluster(reduction, 'rna', 'SLC7A5'), ncol = 4)))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/Treg.pdf'), width = 16, height = 16)
  try(print(wrap_plots(check_marker_cluster(reduction, 'adt', 'CD152'), 
                       check_marker_cluster(reduction, 'adt', 'CD279'), check_marker_cluster(reduction, 'adt', 'CD39'),
                       check_marker_cluster(reduction, 'adt', 'CD103'),
                       check_marker_cluster(reduction, 'adt', 'CD25'), check_marker_cluster(reduction, 'adt', 'CD134'),
                       check_marker_cluster(reduction, 'adt', 'CD137'), check_marker_cluster(reduction, 'adt', 'CD223'), 
                       check_marker_cluster(reduction, 'rna', 'CCR7'), check_marker_cluster(reduction, 'motif', 'MA0850.1', 'FOXP3'), 
                       check_marker_cluster(reduction, 'rna', 'FOXP3'), check_marker_cluster(reduction, 'rna', 'IKZF2'), 
                       check_marker_cluster(reduction, 'rna', 'SMAD3'), check_marker_cluster(reduction, 'rna', 'AHR'), 
                       check_marker_cluster(reduction, 'rna', 'CTLA4'), check_marker_cluster(reduction, 'rna', 'IL2RA'), ncol = 4)))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/Th17.pdf'), width = 20, height = 16)
  try(print(wrap_plots(check_marker_cluster(reduction, 'adt', 'CD161'), check_marker_cluster(reduction, 'adt', 'CD278'),
                       check_marker_cluster(reduction, 'adt', 'CD194'), check_marker_cluster(reduction, 'adt', 'CD196'),
                       check_marker_cluster(reduction, 'rna', 'IL21R'), check_marker_cluster(reduction, 'rna', 'IL23R'), 
                       check_marker_cluster(reduction, 'rna', 'IL17F'), 
                       check_marker_cluster(reduction, 'rna', 'TNF'), check_marker_cluster(reduction, 'rna', 'CCL20'), 
                       check_marker_cluster(reduction, 'motif', 'MA1151.1', 'RORC'),
                       check_marker_cluster(reduction, 'motif', 'MA0071.1', 'RORA'), check_marker_cluster(reduction, 'motif', 'MA0072.1', 'RORA'), check_marker_cluster(reduction, 'rna', 'RORA'), 
                       check_marker_cluster(reduction, 'motif', 'MA1634.1', 'BATF'), check_marker_cluster(reduction, 'rna', 'BATF'), 
                       check_marker_cluster(reduction, 'motif', 'MA1419.1', 'IRF4'), check_marker_cluster(reduction, 'rna', 'IRF4'), 
                       check_marker_cluster(reduction, 'motif', 'MA1419.1', 'IRF4'), check_marker_cluster(reduction, 'rna', 'IRF4'), ncol = 5)))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/Th1.pdf'), width = 16, height = 12)
  try(print(wrap_plots(check_marker_cluster(reduction, 'adt', 'CD183'), check_marker_cluster(reduction, 'adt', 'CD195'),
                       check_marker_cluster(reduction, 'adt', 'CD26'), check_marker_cluster(reduction, 'adt', 'CD94'),
                       check_marker_cluster(reduction, 'adt', 'CD278'), 
                       check_marker_cluster(reduction, 'rna', 'IL18R1'), check_marker_cluster(reduction, 'rna', 'IFNG'), 
                       check_marker_cluster(reduction, 'rna', 'LTA'), check_marker_cluster(reduction, 'rna', 'TNF'), 
                       check_marker_cluster(reduction, 'motif', 'MA0690.1', 'TBX21'), check_marker_cluster(reduction, 'rna', 'TBX21'), ncol = 4)))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/Tfh.pdf'), width = 20, height = 16)
  try(print(wrap_plots(check_marker_cluster(reduction, 'adt', 'CD185'), check_marker_cluster(reduction, 'adt', 'CD84'),
                       check_marker_cluster(reduction, 'adt', 'CD150'), check_marker_cluster(reduction, 'adt', 'CD200'),
                       check_marker_cluster(reduction, 'adt', 'CD272'), check_marker_cluster(reduction, 'adt', 'CD278'), 
                       check_marker_cluster(reduction, 'adt', 'CD279'), check_marker_cluster(reduction, 'adt', 'TIGIT'), 
                       check_marker_cluster(reduction, 'adt', 'CD57'), check_marker_cluster(reduction, 'adt', 'CD154'), 
                       check_marker_cluster(reduction, 'adt', 'CD304'),
                       check_marker_cluster(reduction, 'rna', 'IL6R'), check_marker_cluster(reduction, 'rna', 'TNFSF8'), 
                       check_marker_cluster(reduction, 'rna', 'IL21R'), check_marker_cluster(reduction, 'rna', 'CD40LG'), 
                       check_marker_cluster(reduction, 'motif', 'MA0463.2', 'BCL6'), check_marker_cluster(reduction, 'rna', 'BCL6'), ncol = 5)))
  dev.off()
  
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/Mono.pdf'), width = 16, height = 12)
  try(print(wrap_plots(check_marker_cluster(reduction, 'adt', 'CD14'), check_marker_cluster(reduction, 'adt', 'CD16'),
                       check_marker_cluster(reduction, 'adt', 'CD192'), check_marker_cluster(reduction, 'adt', 'CD83'),
                       check_marker_cluster(reduction, 'adt', 'CD123'), check_marker_cluster(reduction, 'adt', 'CD99'), 
                       check_marker_cluster(reduction, 'adt', 'CD64'), check_marker_cluster(reduction, 'adt', 'CD35'), 
                       check_marker_cluster(reduction, 'adt', 'CD11b'), check_marker_cluster(reduction, 'adt', 'CD88'), 
                       check_marker_cluster(reduction, 'adt', 'CX3CR1'),
                       check_marker_cluster(reduction, 'rna', 'CD83'), ncol = 4)))
  dev.off()
}

p1 <- DimPlot(tcell, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = 'wsnn2_res.0.06')
p2 <- DimPlot(tcell, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")

pdf('../plots/annotation/clusters.pdf', width = 20, height = 8)
print(p1 | p2)
dev.off()

predict_cluster <- function(data, celltype, dataset, cluster){
  reference <- readRDS('../../Time_series/output/reference_qc.RDS')
  reference <- reference[,reference$celltype %in% celltype]
  
  DefaultAssay(reference) <- "RNA"
  reference <- SCTransform(reference, method = "glmGamPoi", vars.to.regress = "percent.mt") %>% RunPCA(reduction.name = "pca") %>%
    RunHarmony(group.by.vars = c("sample"), max.iter.harmony = 40, reduction = 'pca', assay.use = 'SCT',project.dim = FALSE,  reduction.save = "harmony_RNA")
  reference <- RunUMAP(reference, reduction = "harmony_RNA", dims = 1:30, assay = 'SCT', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_', return.model = T)
  
  DefaultAssay(reference) <- 'ADT'
  # we will use all ADT features for dimensional reduction
  # we set a dimensional reduction name to avoid overwriting the
  VariableFeatures(reference) <- rownames(reference[["ADT"]])
  reference <- NormalizeData(reference, normalization.method = 'CLR', margin = 2) %>%
    ScaleData() %>% RunPCA(reduction.name = 'apca') %>%
    RunHarmony(group.by.vars = c("sample", "cryopreservation"), max.iter.harmony = 40, reduction = 'apca', assay.use = 'ADT',project.dim = FALSE,  reduction.save = "harmony_ADT")
  reference <- RunUMAP(reference, reduction = "harmony_ADT", dims = 1:30, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_', return.model = T)
  
  DefaultAssay(reference) <- "peaks"
  reference <- RunTFIDF(reference)
  reference <- FindTopFeatures(reference, min.cutoff = 'q0')
  reference <- RunSVD(reference)
  reference <- RunHarmony(reference, group.by.vars = c("sample"), max.iter.harmony = 40, reduction = 'lsi', assay.use = 'peaks',project.dim = FALSE,  reduction.save = "harmony_peaks")
  reference <- RunUMAP(reference, reduction = "harmony_peaks", dims = 2:30, assay = 'peaks', reduction.name = 'atac.umap', reduction.key = 'atacUMAP_', return.model = T)
  
  reference <- FindMultiModalNeighbors(object = reference,
                                  reduction.list = list("harmony_RNA", "harmony_peaks", "harmony_ADT"),
                                  dims.list = list(1:30, 2:30, 1:30))
  reference <- RunUMAP(reference, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", return.model = T)
  
  reference <- RunSPCA(reference, assay = 'SCT', graph = 'wsnn', reduction.name = 'rna.spca', reduction.key = 'rnaSPC')
  reference <- RunSPCA(reference, assay = 'ADT', graph = 'wsnn', reduction.name = 'adt.spca', reduction.key = 'adtSPC')
  
  DefaultAssay(data) <- 'SCT'
  DefaultAssay(reference) <- 'SCT'
  
  anchors.rna <- FindTransferAnchors(
    reference = reference,
    query = data,
    normalization.method = "SCT",
    reference.reduction = "rna.spca",
    reduction = 'pcaproject',
    dims = 1:30,
    reference.assay = 'SCT',
    query.assay = 'SCT'
  )
  
  anchors.adt <- FindTransferAnchors(
    reference = reference,
    query = data,
    reference.reduction = "adt.spca",
    reduction = 'pcaproject',
    dims = 1:30,
    reference.assay = 'ADT',
    query.assay = 'ADT'
  )
  
  anchors.atac <- FindTransferAnchors(
    reference = reference,
    query = data,
    reference.reduction = "lsi",
    reduction = 'lsiproject',
    dims = 2:30,
    reference.assay = 'peaks',
    query.assay = 'peaks'
  )
  
  gc()
  
  predictions.rna <- TransferData(anchorset = anchors.rna, refdata = reference$celltype, weight.reduction = "pcaproject")
  predictions.adt <- TransferData(anchorset = anchors.adt, refdata = reference$celltype, weight.reduction = "pcaproject")
  predictions.atac <- TransferData(anchorset = anchors.atac, refdata = reference$celltype, weight.reduction = "lsiproject")
  
  gc()
  
  predictions.max <- data.table(rna = rowMax(as.matrix(predictions.rna[,-1])),
                                adt = rowMax(as.matrix(predictions.adt[,-1])),
                                atac = rowMax(as.matrix(predictions.atac[,-1])))
  rownames(predictions.max) <- rownames(predictions.rna)
  predictions.max.label <- apply(predictions.max, 1, which.max)
  names(predictions.max.label) <- rownames(predictions.rna)
  
  predictions <- rbind(predictions.rna[names(predictions.max.label[predictions.max.label == 1]),],
                       predictions.adt[names(predictions.max.label[predictions.max.label == 2]),],
                       predictions.atac[names(predictions.max.label[predictions.max.label == 3]),])
  predictions <- predictions[rownames(predictions.rna),1]
  
  data$annotated_predicted.celltype <- factor(predictions, levels = celltype)
  
  data$annotated_predicted.celltype.rna <- factor(predictions.rna$predicted.id, levels = celltype)
  
  data$annotated_predicted.celltype.adt <- factor(predictions.adt$predicted.id, levels = celltype)
  
  data$annotated_predicted.celltype.atac <- factor(predictions.atac$predicted.id, levels = celltype)
  
  p1 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = 'wnn_clusters')
  p2 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/predicted_annotated_predicted_WNN.pdf'), width = 20, height = 8)
  print(p1 | p2)
  dev.off()
  
  p1 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = 'rna_clusters')
  p2 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/predicted_annotated_predicted_RNA.pdf'), width = 20, height = 8)
  print(p1 | p2)
  dev.off()
  
  p1 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = 'adt_clusters')
  p2 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/predicted_annotated_predicted_ADT.pdf'), width = 20, height = 8)
  print(p1 | p2)
  dev.off()
  
  p1 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = 'atac_clusters')
  p2 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/predicted_annotated_predicted_ATAC.pdf'), width = 20, height = 8)
  print(p1 | p2)
  dev.off()
  
  p1 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = 'wnn2_clusters')
  p2 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = "annotated_predicted.celltype")
  pdf(paste0('../plots/', dataset, '/tcell/', cluster, '/predicted_annotated_predicted_WNN2.pdf'), width = 20, height = 8)
  print(p1 | p2)
  dev.off()
  
  saveRDS(data, file = paste0('../output/', dataset, '/predicted_tcell_', cluster, '.RDS'))
}

#CD4+ naive
analysis(tcell[, tcell$wsnn2_res.0.06 %in% c(0)], 'annotation', 'cluster0')

data <- readRDS('../output/annotation/tcell_cluster0.RDS')
predict_cluster(data, c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)'), 'annotation', 'cluster0')

data <- readRDS('../output/annotation/predicted_tcell_cluster0.RDS')

#CD8+ naive
analysis(tcell[, tcell$wsnn2_res.0.06 %in% c(2)], 'annotation', 'cluster2')

data <- readRDS('../output/annotation/tcell_cluster2.RDS')
predict_cluster(data, c('CD8+ Naive (Resting)', 'CD8+ Naive (Activated)',
                        'CD8+ Regulatory'), 'annotation', 'cluster2')

data <- readRDS('../output/annotation/predicted_tcell_cluster2.RDS')

#CD8+ memory
analysis(tcell[, tcell$wsnn2_res.0.06 %in% c(3,4,6)], 'annotation', 'cluster3')

data <- readRDS('../output/annotation/tcell_cluster3.RDS')
predict_cluster(data, c('CD8+ Memory (Resting)', 'CD8+ Memory (Activated)',
                        'MAITs (Resting)', 'MAITs (Activated)', 'Gamma Delta'), 'annotation', 'cluster3')

data <- readRDS('../output/annotation/predicted_tcell_cluster3.RDS')

#CD4+ memory
analysis(tcell[, tcell$wsnn2_res.0.06 %in% c(1)], 'annotation', 'cluster1')

data <- readRDS('../output/annotation/tcell_cluster1.RDS')
predict_cluster(data, c('CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
                        'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1', 
                        'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17', 
                        'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other'), 'annotation', 'cluster1')

data <- readRDS('../output/annotation/predicted_tcell_cluster1.RDS')

#Tfh
analysis(tcell[, tcell$wsnn2_res.0.06 %in% c(5)], 'annotation', 'cluster5')

data <- readRDS('../output/annotation/tcell_cluster5.RDS')
predict_cluster(data, c('CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh'), 'annotation', 'cluster5')

data <- readRDS('../output/annotation/predicted_tcell_cluster5.RDS')

cluster0 <- readRDS('../output/annotation/predicted_tcell_cluster0.RDS')
cluster1 <- readRDS('../output/annotation/predicted_tcell_cluster1.RDS')
cluster2 <- readRDS('../output/annotation/predicted_tcell_cluster2.RDS')
cluster3 <- readRDS('../output/annotation/predicted_tcell_cluster3.RDS')
cluster5 <- readRDS('../output/annotation/predicted_tcell_cluster5.RDS')

tcell$celltype <- NA
tcell$celltype[colnames(cluster0)] <- as.character(cluster0$annotated_predicted.celltype)
tcell$celltype[colnames(cluster1)] <- as.character(cluster1$annotated_predicted.celltype)
tcell$celltype[colnames(cluster2)] <- as.character(cluster2$annotated_predicted.celltype)
tcell$celltype[colnames(cluster3)] <- as.character(cluster3$annotated_predicted.celltype)
tcell$celltype[colnames(cluster5)] <- as.character(cluster5$annotated_predicted.celltype)

DimPlot(tcell, reduction = 'wnn2.umap', group.by = 'celltype')

tcell$celltype <- factor(tcell$celltype, levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
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

data <- tcell
rm(tcell)

Idents(data) <- 'celltype'
data <- RunUMAP(data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", return.model = T)

data <- RunSPCA(data, assay = 'SCT', graph = 'wsnn', reduction.name = 'rna.spca', reduction.key = 'rnaSPC')
data <- RunSPCA(data, assay = 'ADT', graph = 'wsnn', reduction.name = 'adt.spca', reduction.key = 'adtSPC')

saveRDS(data, '../output/tcell_downsample_annotated.RDS')

p1 <- DimPlot(data, group.by = 'celltype', reduction = 'wnn2.umap') + guides(col = guide_legend(ncol = 1, override.aes = list(size=4)))

pdf('../plots/annotation/annotated.pdf', width = 10, height = 8)
p1
dev.off()

p1 <- DimPlot(data, group.by = 'celltype', reduction = 'wnn.umap') + guides(col = guide_legend(ncol = 1, override.aes = list(size=4)))
p2 <- DimPlot(data, group.by = 'celltype', reduction = 'rna.umap') + guides(col = guide_legend(ncol = 1, override.aes = list(size=4)))
p3 <- DimPlot(data, group.by = 'celltype', reduction = 'adt.umap') + guides(col = guide_legend(ncol = 1, override.aes = list(size=4)))
p4 <- DimPlot(data, group.by = 'celltype', reduction = 'atac.umap') + guides(col = guide_legend(ncol = 1, override.aes = list(size=4)))
p5 <- DimPlot(data, group.by = 'celltype', reduction = 'wnn2.umap') + guides(col = guide_legend(ncol = 1, override.aes = list(size=4)))

pdf('../plots/annotation/annotated_all.pdf', width = 30, height = 16)
wrap_plots(p1, p2, p3, p4, p5, ncol = 3)
dev.off()
