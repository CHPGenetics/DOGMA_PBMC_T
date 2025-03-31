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

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')

load('../output/data_harmony_qc_tcell_updated.RData')

plan('multicore')
options(future.globals.maxSize= 60*1024^3)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

DefaultAssay(data) <- "peaks"
peaks <- CallPeaks(data, macs2.path = "/ihome/wchen/zhongli/.conda/envs/MACS2/bin/macs2")
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

save(peaks, file = '../data/tcell_peaks.RData')

macs2_counts <- FeatureMatrix(
  fragments = Fragments(data),
  features = peaks,
  cells = colnames(data)
)

save(macs2_counts, file = '../data/tcell_macs2_counts.RData')

# create a new assay using the MACS2 peak set and add it to the Seurat object
data[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(data),
  annotation = annotations
)

save(data, file = '../output/tcell.RData')

DefaultAssay(data) <- "peaks"
data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = 'q0')
data <- RunSVD(data)
data <- RunHarmony(data, group.by.vars = c("sample"), max.iter.harmony = 20, reduction = 'lsi', assay.use = 'peaks',project.dim = FALSE,  reduction.save = "harmony_peaks")
data <- RunUMAP(data, reduction = "harmony_peaks", dims = 2:30, assay = 'peaks', reduction.name = 'atac.umap', reduction.key = 'atacUMAP_')
data <- FindNeighbors(data, dims = 2:30, reduction = 'harmony_peaks')
data <- FindClusters(data, resolution = 0.2, graph.name = 'peaks_snn', algorithm = 3)
data$harmony_atac_clusters <- Idents(data)

# Now run multimodal neighbors and embedding
data <- FindMultiModalNeighbors(object = data,
                                reduction.list = list("harmony_RNA", "harmony_peaks", "harmony_ADT"),
                                dims.list = list(1:30, 2:30, 1:30))
data <- RunUMAP(data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_" )
data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, resolution = 0.2)

data$harmony_wnn_clusters <- Idents(data)

data <- FindMultiModalNeighbors(object = data,
                                reduction.list = list("harmony_peaks", "harmony_ADT"),
                                dims.list = list(2:30, 1:30),
                                knn.graph.name = "wknn2",
                                snn.graph.name = "wsnn2",
                                weighted.nn.name = "weighted.nn2")
data <- RunUMAP(data, nn.name = "weighted.nn2", reduction.name = "wnn2.umap", reduction.key = "wnn2UMAP_" )
data <- FindClusters(data, graph.name = "wsnn2", algorithm = 3, resolution = 0.2)
data$harmony_wnn2_clusters <- Idents(data)

save(data, file = '../output/tcell.RData')

dataset <- 'tcell'

p1 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = 'harmony_wnn_clusters')
p2 <- DimPlot(data, reduction = 'wnn.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1")
pdf(paste0('../plots/', dataset, '/predicted_WNN.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = 'harmony_rna_clusters')
p2 <- DimPlot(data, reduction = 'rna.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1")
pdf(paste0('../plots/', dataset, '/predicted_RNA.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = 'harmony_adt_clusters')
p2 <- DimPlot(data, reduction = 'adt.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1")
pdf(paste0('../plots/', dataset, '/predicted_ADT.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = 'harmony_atac_clusters')
p2 <- DimPlot(data, reduction = 'atac.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1")
pdf(paste0('../plots/', dataset, '/predicted_ATAC.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

p1 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = 'harmony_wnn2_clusters')
p2 <- DimPlot(data, reduction = 'wnn2.umap', label = TRUE, repel = TRUE, group.by = "predicted.celltype.l1")
pdf(paste0('../plots/', dataset, '/predicted_WNN2.pdf'), width = 20, height = 8)
print(p1 | p2)
dev.off()

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(data) <- "peaks"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(data), pwm = pwm_set, genome = 'BSgenome.Hsapiens.UCSC.hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
data <- SetAssayData(data, assay = 'peaks', slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes
data <- RunChromVAR(
  object = data,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

save(data, file = '../output/tcell_chrom.RData')

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
p1 <- FeaturePlot(data, 'Activated1', reduction = 'wnn.umap', min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle('Activated')
p2 <- FeaturePlot(data, 'Resting1', reduction = 'wnn.umap', min.cutoff = 'q1', max.cutoff = 'q99') + ggtitle('Resting') 
pdf('../plots/tcell/activated_resting.pdf', width = 16, height = 8)
print(p1 | p2)
dev.off()

save(data, file = '../output/tcell_chrom.RData')