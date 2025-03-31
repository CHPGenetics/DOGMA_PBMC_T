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
library(BiocFileCache)
library(SummarizedExperiment)
library(cowplot)
library(ggrepel)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

data <- readRDS('../output/tcell_annotated_updated.RDS')

Idents(data) <- 'celltype_updated'

motif <- readRDS('../output/ArchR/DOGMA_filtered_multiome/motif.RDS')
motif <- assays(motif)[[2]]
colnames(motif) <- gsub('#', '_', colnames(motif))
motif <- motif[,colnames(data)]
rownames(motif) <- gsub('_', '-', rownames(motif))
data[['chromvar2']] <- CreateAssayObject(data = motif)
rm(motif)

rna <- readRDS('../output/ArchR/DOGMA_filtered_multiome/rna.RDS')
rowdata <- rowData(rna)
rna <- assays(rna)[[1]]
colnames(rna) <- gsub('#', '_', colnames(rna))
rna <- rna[,colnames(data)]
rownames(rna) <- rowdata$name
data[['RNA2']] <- CreateAssayObject(data = rna)
rm(rna, rowdata)

corGEM_MM <- read.csv('../output/ArchR/DOGMA_filtered_multiome/corGEM_MM.csv', row.names = 'X')
regulator <- corGEM_MM[corGEM_MM$TFRegulator != 'No',]
motif <- data@assays$chromvar2@data[gsub('_', '-', regulator$MotifMatrix_name),]
data[['chromvar2_sub']] <- CreateAssayObject(data = motif)

motif.markers <- FindAllMarkers(data, assay = 'chromvar2_sub', logfc.threshold = 0, min.pct = 0, test.use = 'LR')
write.xlsx(motif.markers, file = '../plots/regulator/heatmap_motif_regulator.xlsx', colNames = T, rowNames = T, overwrite = T)
motif.top10 <- motif.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf('../plots/regulator/heatmap_motif_regulator.pdf', width = 16, height = 16)
print(DoHeatmap(data, assay = 'chromvar2_sub', slot = 'data', features = motif.top10$gene, angle = 90, size = 4) + NoLegend())
dev.off()

rna <- data@assays$RNA2@data[regulator$GeneExpressionMatrix_name,]
rownames(motif) <- str_split_fixed(rownames(motif), '-', 2)[,1]
rownames(rna) <- str_split_fixed(rownames(rna), '-', 2)[,1]
table(rownames(motif) == rownames(rna))
# rna <- t(scale(t(rna)))

data_plot <- as.data.frame(matrix(rep(NA, 4), ncol = 4))
colnames(data_plot) <- c('rna', 'motif', 'celltype', 'name')
for(i in levels(data$celltype_updated)){
  barcodes <- colnames(data)[data$celltype_updated == i]
  df <- data.frame(rna = rowMeans(rna[,barcodes]),
                   motif = rowMeans(motif[,barcodes]),
                   celltype = i,
                   name = rownames(rna))
  data_plot <- rbind(data_plot, df)
}
data_plot <- data_plot[-1,]
data_plot$celltype <- factor(data_plot$celltype, levels = levels(data$celltype_updated))

data_plot_ht <- as.data.frame(matrix(rep(NA, 52), ncol = 1))
for(i in levels(data$celltype_updated)){
  barcodes <- colnames(data)[data$celltype_updated == i]
  data_plot_ht <- cbind(data_plot_ht, rowMeans(motif[,barcodes]))
}
data_plot_ht <- data_plot_ht[,-1]
colnames(data_plot_ht) <- levels(data$celltype_updated)
data_plot_ht <- as.matrix(data_plot_ht)

library(pheatmap)
pdf('../plots/regulator/heatmap_motif_regulator_order.pdf', width = 16, height = 16)
p <- pheatmap(data_plot_ht, cluster_rows = T, cluster_cols = F)
p
dev.off()

order <- p$tree_row$labels[p$tree_row$order]
pheatmap(data_plot_ht[order[c(1:10, 14:20, 33:48)],], cluster_rows = F, cluster_cols = F)
p1 <- pheatmap(data_plot_ht[order[-c(1:10, 14:20, 33:48)],], cluster_rows = T, cluster_cols = F)
order1 <- p1$tree_row$labels[p1$tree_row$order]
pheatmap(data_plot_ht[c(order[c(1:10, 14:20, 33:48)], order1),], cluster_rows = F, cluster_cols = F)

gene_order <- c(order[c(1:10, 14:20, 33:48)], order1)

# gene_order <- motif.top10$gene
# gene_order <- str_split_fixed(gene_order, '-', 2)[,1]
# gene_order <- gene_order[!duplicated(gene_order)]
data_plot$name <- factor(data_plot$name, levels = gene_order)

data_plot_new <- as.data.frame(matrix(rep(NA, 5), ncol = 5))
colnames(data_plot_new) <- c('rna', 'motif', 'celltype', 'name', 'rna_scaled')
for(i in names(table(data_plot$name))){
  df <- data_plot[data_plot$name == i,]
  df$rna_scaled <- scale(df$rna)
  data_plot_new <- rbind(data_plot_new, df)
}
data_plot_new <- data_plot_new[-1,]
data_plot_new$celltype <- factor(data_plot_new$celltype, levels = levels(data$celltype_updated))
data_plot_new$name <- factor(data_plot_new$name, levels = gene_order)

pdf('../plots/regulator/motif_regulator.pdf', width = 14, height = 7)
ggplot(data_plot_new, aes(y = celltype, x = name, size = rna_scaled, fill = motif))+
  geom_point(shape = 21)+
  scale_size(range = c(0,7))+
  scale_fill_gradient2(low = 'dodgerblue3', mid = 'white', high = 'firebrick')+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'bottom')+
  labs(x = '', y = '', fill = 'Motif deviation', size = 'TF expression') + 
  guides(size = guide_legend(nrow = 1))
dev.off()

# data_plot_hc <- as.data.frame(matrix(NA, nrow = 51, ncol = 1))
# colnames(data_plot_hc) <- levels(data$celltype)
# for(i in levels(data$celltype)){
#   barcodes <- colnames(data)[data$celltype == i]
#   df <- rowMeans(motif[,barcodes])
#   data_plot_hc <- cbind(data_plot_hc, df)
# }
# data_plot_hc <- data_plot_hc[,-1]
# colnames(data_plot_hc) <- levels(data$celltype)
# data_plot_hc <- as.matrix(data_plot_hc)
# 
# pheatmap(data_plot_hc, cluster_rows = T, cluster_cols = F)

corGEM_MM$name <- str_split_fixed(corGEM_MM$GeneExpressionMatrix_name, '-', 2)[,1]
p <- ggplot(data.frame(corGEM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() +
  geom_text_repel(color = 'black', segment.color = 'black', min.segment.length = 0, max.overlaps = Inf, aes(label=ifelse(name %in% c('LEF1', 'CEBPD', 'RORC', 'ETS1', 'ELF2', 'IKZF1', 'NFKB1', 'REL', 'MAF'),name,'')))+
  theme_cowplot() +
  geom_vline(xintercept = 0, lty = "dashed") +
  scale_color_manual(values = c("Negative" = "dodgerblue3", "No"="darkgrey", "Positive"="firebrick3")) +
  xlab("Correlation with TF expression") +
  ylab("Max TF motif delta") +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, max(corGEM_MM$maxDelta)*1.05)
  )+labs(color = 'TF regulators')+theme(legend.position = 'bottom')

pdf('../plots/regulator/Correlation-Motif-Deviation-TF-Expression.pdf', width = 5, height = 5)
p
dev.off()

p <- ggplot(data.frame(corGEM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() +
  geom_text_repel(color = 'black', segment.color = 'black', min.segment.length = 0, max.overlaps = Inf, aes(label=ifelse(name %in% c('LEF1', 'CEBPD', 'RORC', 'ETS1', 'ELF2', 'IKZF1', 'NFKB1', 'REL', 'JUND', 'TCF7L2', 'ETS2', 'CREB1'),name,'')))+
  theme_cowplot() +
  geom_vline(xintercept = 0, lty = "dashed") +
  scale_color_manual(values = c("Negative" = "dodgerblue3", "No"="darkgrey", "Positive"="firebrick3")) +
  xlab("Correlation with TF expression") +
  ylab("Max TF motif delta") +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, max(corGEM_MM$maxDelta)*1.05)
  )+labs(color = 'TF regulators')+theme(legend.position = 'bottom')

pdf('../plots/regulator/Correlation-Motif-Deviation-TF-Expression-updated.pdf', width = 5, height = 5)
p
dev.off()

df = fread('../output/ArchR/DOGMA_filtered_multiome/Plots/Plot-MonocleCD4Naive_fine-Traj-Heatmap_Cor.csv')
table(regulator$GeneExpressionMatrix_name %in% df$x)
