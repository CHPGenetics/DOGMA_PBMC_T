library(Signac)
library(Seurat)
library(SummarizedExperiment)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

data <- readRDS('../output/tcell_annotated_updated.RDS')
load('../output/ArchR/DOGMA_filtered_multiome/PeakMatrix.RData')

peaks_mat_df <- data.frame(chr=seqnames(peaks_mat),
                           start=start(peaks_mat),
                           end=end(peaks_mat),
                           distToTSS = rowRanges(peaks_mat)$distToTSS,
                           GC = rowRanges(peaks_mat)$GC)
rownames(peaks_mat_df) <- GRangesToString(rowRanges(peaks_mat), sep = c("-", "-"))

meta <- as.data.frame(colData(peaks_mat))
mat <- assay(peaks_mat)

rownames(mat) <- rownames(peaks_mat_df)
colnames(mat) <- gsub('#', '_', colnames(mat))

DefaultAssay(data) <- 'peaks'
data[['peaks2']] <- CreateChromatinAssay(counts = mat,
                                         fragments = Fragments(data),
                                         annotation = Annotation(data))
rm(mat, meta, peaks_mat, peaks_mat_df)

DefaultAssay(data) <- 'peaks2'
peak_annotation <- ClosestFeature(data, regions = rownames(data))
rownames(peak_annotation) <- peak_annotation$query_region

rna <- data@assays$RNA@counts
adt <- data@assays$ADT@counts
peak <- data@assays$peaks2@counts

saveRDS(peak_annotation, file = '../output/pseudo_bulk/peak_annotation.RDS')
saveRDS(rna, file = '../output/pseudo_bulk/rna.RDS')
saveRDS(adt, file = '../output/pseudo_bulk/adt.RDS')
saveRDS(peak, file = '../output/pseudo_bulk/peak.RDS')
