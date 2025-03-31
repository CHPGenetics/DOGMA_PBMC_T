rm(list=ls())
gc()

# Load Packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(biomaRt)
library(httr)
library(rtracklayer)
library(liftOver)
library(GenomicRanges)

# Parameters
DATA_PATH <- "/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Input/"

# Load Data
tcell_dogma <- readRDS(paste0(DATA_PATH, "tcell_annotated_updated.RDS"))

# Gene
ensembl_list <- rownames(tcell_dogma@assays$RNA@counts)
set_config(config(ssl_verifypeer = 0L))
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 37)
gene_bp <- getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"), 
                 filters = "external_gene_name", values = ensembl_list, mart = ensembl)
gene_bp_filter <- gene_bp[which(gene_bp$chromosome_name %in% as.character(1:22)),]
gene_bp_filter$width <- gene_bp_filter$end_position - gene_bp_filter$start_position
gene_bp_filter <- arrange(gene_bp_filter, desc(width))
gene_bp_filter <- gene_bp_filter[!duplicated(gene_bp_filter$external_gene_name), ]
genes <- gene_bp_filter[,c("external_gene_name","chromosome_name", "start_position","end_position")]
genes <- genes %>% arrange(chromosome_name, start_position)
colnames(genes) <- c("GENE", "CHR", "START", "END")
write.table(genes$GENE, row.names = F, col.names = F, quote = F, sep = "\t",
            file = paste0(WORK_PATH, "Gene_list.txt"))
write.table(genes, row.names = F, col.names = T, quote = F, sep = "\t",
            file = paste0(WORK_PATH, "Gene_coord_all.txt"))

## DE
DefaultAssay(tcell_dogma) <- "SCT"
n_top <- round(nrow(genes)/10)
all_clusters <- unique(tcell_dogma$celltype_updated)
gene_set <- list()
for (c in 1:length(all_clusters)) {
  ident_x <- all_clusters[c]
  ident_y <- setdiff(all_clusters, ident_x)
  gene_set[[c]] <- FindMarkers(tcell_dogma, 
                               ident.1 = ident_x,
                               ident.2 = ident_y,
                               logfc.threshold = 0,
                               min.pct = 0.001,
                               only.pos = T)
  gene_set_f <- gene_set[[c]] %>% 
    top_n(n = -n_top, wt = p_val)
  gene_list <- rownames(gene_set_f)
  gene_list <- intersect(gene_list, genes$GENE)
  write.table(gene_list,
              file = paste0(WORK_PATH, "Gene_list_", gsub(" ", "_", all_clusters[c]), ".txt"),
              sep = "\t", quote = F, col.names =F, row.names =F)
}

# ATAC
## Get peaks for hg38
list <- rownames(tcell_dogma@assays$peaks@data)
peaks <- data.frame(GENE = list)
peaks <- peaks %>%
  mutate(GENE_copy = GENE) %>% 
  separate(col = GENE_copy, into = c("CHR", "START", "END"), sep = "-") %>%
  mutate(CHR = gsub("chr", "", CHR)) %>%
  dplyr::filter(CHR %in% 1:22) %>%
  arrange(CHR, START)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
chain = import.chain(path)
gr <- StringToGRanges(regions = peaks$GENE, sep = c("-", "-"))
gr.hg19 <- liftOver(gr, chain)
names(gr.hg19) <- peaks$GENE
correspondence <- elementNROWS(gr.hg19)
gr.hg19 <- gr.hg19[correspondence == 1]
gr.hg19 <- unlist(gr.hg19) %>% as.data.frame()
gr.hg19$GENE <- rownames(gr.hg19)
peaks <- inner_join(peaks, gr.hg19, by = "GENE")
peaks <- data.frame(GENE = peaks$GENE,
                    CHR = peaks$CHR,
                    START = peaks$start,
                    END = peaks$end)
write.table(peaks$GENE, row.names = F, col.names = F, quote = F, sep = "\t",
            file = paste0(WORK_PATH, "Peak_list.txt"))
write.table(peaks, row.names = F, col.names = T, quote = F, sep = "\t",
            file = paste0(WORK_PATH, "Peak_coord_all.txt"))

## DA
DefaultAssay(tcell_dogma) <- "peaks"
n_top <- round(nrow(peaks)/10)
all_clusters <- unique(tcell_dogma$celltype_updated)
peak_set <- list()
for (c in 1:length(all_clusters)) {
  ident_x <- all_clusters[c]
  ident_y <- setdiff(all_clusters, ident_x)
  peak_set[[c]] <- FindMarkers(tcell_dogma, 
                               ident.1 = ident_x,
                               ident.2 = ident_y,
                               logfc.threshold = 0,
                               min.pct = 0.001,
                               only.pos = T)
  peak_set_f <- peak_set[[c]] %>% top_n(n = -n_top, wt = p_val)
  peak_list <- rownames(peak_set_f)
  peak_list <- intersect(peak_list, peaks$GENE)
  write.table(peak_list,
              file = paste0(WORK_PATH, "Peak_list_", gsub(" ", "_", all_clusters[c]), ".txt"),
              sep = "\t", quote = F, col.names =F, row.names =F)
}

## Peak2
peak2 <- readRDS('/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/pseudo_bulk/peak.RDS')
list <- rownames(peak2)
rm(peak2)
peaks <- data.frame(GENE = list)
peaks <- peaks %>%
  mutate(GENE_copy = GENE) %>% 
  separate(col = GENE_copy, into = c("CHR", "START", "END"), sep = "-") %>%
  mutate(CHR = gsub("chr", "", CHR)) %>%
  dplyr::filter(CHR %in% 1:22) %>%
  arrange(CHR, START)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
chain = import.chain(path)
gr <- StringToGRanges(regions = peaks$GENE, sep = c("-", "-"))
gr.hg19 <- liftOver(gr, chain)
names(gr.hg19) <- peaks$GENE
correspondence <- elementNROWS(gr.hg19)
gr.hg19 <- gr.hg19[correspondence == 1]
gr.hg19 <- unlist(gr.hg19) %>% as.data.frame()
gr.hg19$GENE <- rownames(gr.hg19)
peaks <- inner_join(peaks, gr.hg19, by = "GENE")
peaks <- data.frame(GENE = peaks$GENE,
                    CHR = peaks$CHR,
                    START = peaks$start,
                    END = peaks$end)
write.table(peaks$GENE, row.names = F, col.names = F, quote = F, sep = "\t",
            file = paste0(WORK_PATH, "Peak2_list.txt"))
write.table(peaks, row.names = F, col.names = T, quote = F, sep = "\t",
            file = paste0(WORK_PATH, "Peak2_coord_all.txt"))






