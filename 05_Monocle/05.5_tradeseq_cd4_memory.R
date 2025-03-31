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
library(monocle3)
library(SeuratWrappers)
library(magrittr)
library(pheatmap)
library(SummarizedExperiment)
library(tradeSeq)
library(slingshot)
library(condiments)
library(tidyr)
library(cowplot)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

seurat <- readRDS('../output/monocle_2/tcell_annotated_updated_2_conditions_cd4_memory.RDS')

seurat@reductions$umap <- seurat@reductions$atac.umap
seurat@reductions$umap@key <- 'UMAP_'
# 
# data <- as.SingleCellExperiment(seurat, assay = "peaks")
# 
# data$condition <- factor(as.character(data$condition), levels = c("Act_IL1B_IL23","Act_IL1B_IL23_PGE2"))
# data$celltype <- factor(as.character(data$celltype_updated), levels = levels(data$celltype_updated)[5:12])
# theme_set(theme_classic())
# 
# df <- bind_cols(
#   as.data.frame(reducedDims(data)$UMAP),
#   as.data.frame(colData(data))
# )
# p1 <- ggplot(df, aes(x = atacUMAP_1, y = atacUMAP_2, col = condition)) +
#   geom_point(size = .7) +
#   labs(col = "condition")
# p1
# 
# p2 <- ggplot(df, aes(x = atacUMAP_1, y = atacUMAP_2, col = celltype)) +
#   geom_point(size = .7) +
#   labs(col = "celltype")
# p2
# 
# data <- condiments::imbalance_score(data, conditions = data$condition,
#                                     k = 20, smooth = 40, dimred = "UMAP")
# df$scores <- data$scores$scaled_scores
# ggplot(df, aes(x = atacUMAP_1, y = atacUMAP_2, col = scores)) +
#   geom_point(size = .7) +
#   scale_color_viridis_c(option = "C") +
#   labs(col = "Scores")
# 
# library(slingshot)
# data <- slingshot(data, reducedDim = 'UMAP', clusterLabels = data$celltype,
#                   start.clus = 'CD4+ Memory (Resting) - Other', approx_points = 150)
# # data@int_metadata$slingshot <- data$slingshot
# mst <- slingMST(data, as.df = TRUE)
# p3 <- ggplot(df, aes(x = atacUMAP_1, y = atacUMAP_2)) +
#   geom_point(size = 1, alpha = .65, col = "grey70", shape = 21,
#              aes(fill = celltype)) +
#   labs(fill = "Cell Clusters") +
#   geom_point(size = 3, data = mst) +
#   geom_path(size = 1.5, data = mst, aes(group = Lineage))
# p3
# 
# set.seed(100)
# BiocParallel::bpparam()
# topologyTest(data, conditions = data$condition, rep = 10, parallel = T)
# 
# sdss <- slingshot_conditions(data, data$condition, approx_points = FALSE,
#                              extend = "n", reweight = FALSE, reassign = FALSE)
# 
# msts <- lapply(sdss, slingMST, as.df = TRUE) %>%
#   bind_rows(.id = "condition") %>%
#   arrange(condition)
# 
# p4 <- ggplot(df, aes(x = atacUMAP_1, y = atacUMAP_2, col = condition)) +
#   geom_point(size = .7, alpha = .1) +
#   scale_color_brewer(palette = "Accent") +
#   geom_point(data = msts, size = 3) +
#   geom_path(data = msts, aes(group = interaction(Lineage, condition)),
#             size = 2)
# p4
# 
# lineages <- lapply(sdss, slingCurves, as.df = TRUE) %>%
#   bind_rows(.id = "condition") %>%
#   arrange(Order)
# 
# position <- data.frame(
#   "atacUMAP_1" = c(-1, -4, -3, -3.5),
#   "atacUMAP_2" = c(-5, -1, -5, -0.5),
#   "condition" = "Act_IL1b_IL23",
#   "text" = paste0("Lineage ", 1:4)
# )
# 
# p4_1 <- ggplot(df, aes(x = atacUMAP_1, y = atacUMAP_2, col = condition)) +
#   geom_point(size = .7, alpha = .2) +
#   scale_color_brewer(palette = "Accent") +
#   geom_path(data = lineages, size = 1.5, aes(group = interaction(Lineage, condition)))+
#   geom_text(data = position, col = "black", aes(label = text), size = 5)
# p4_1
# 
# mapping <- matrix(rep(1:4, each = 2), nrow = 4, ncol = 2, byrow = TRUE)
# mapping
# 
# sds <- merge_sds(sdss[[1]], sdss[[2]],
#                  condition_id = names(sdss), mapping = mapping)
# 
# df <- cbind(df, slingPseudotime(sds)[rownames(df),])
# 
# ggplot(df, aes(x = atacUMAP_1, y = atacUMAP_2, col = Lineage1)) +
#   geom_point(size = .7) +
#   scale_color_viridis_c() +
#   labs(col = "Pseudotime")
# ggplot(df, aes(x = atacUMAP_1, y = atacUMAP_2, col = Lineage3)) +
#   geom_point(size = .7) +
#   scale_color_viridis_c() +
#   labs(col = "Pseudotime")
# 
# df$cells <- rownames(df)
# df_new <- full_join(
#   df %>%
#     dplyr::select(cells, atacUMAP_1, atacUMAP_2, celltype, condition),
#   slingPseudotime(sds) %>%
#     as.data.frame() %>%
#     mutate(cells = rownames(.))
# ) %>%
#   pivot_longer(starts_with("Lineage"), names_to = "Curve", values_to = "pst")
# 
# p5 <- ggplot(df_new, aes(x = pst)) +
#   geom_density(alpha = .4, aes(fill = condition), col = "transparent") +
#   geom_density(aes(col = condition), fill = "transparent", size = 1.5) +
#   guides(col = FALSE) +
#   scale_fill_brewer(palette = "Accent") +
#   scale_color_brewer(palette = "Accent") +
#   labs(x = "Pseudotime", fill = "Type") +
#   facet_wrap(~ Curve, scales = "free_x")
# p5
# 
# progressionTest(sds, conditions = data[,cellnames(sds)]$condition, lineages = TRUE)
# 
# differentiationTest(sds, conditions = data[,cellnames(sds)]$condition, pairwise = TRUE)
# 
# weights <- condiments:::.sling_reassign(sds)
# df_new <- df_new %>%
#   full_join(weights %>%
#               as.data.frame() %>%
#               mutate(cells = rownames(.)) %>%
#               dplyr::rename("Lineage1" = V1, "Lineage2" = V2, "Lineage3" = V3, "Lineage4" = V4) %>%
#               pivot_longer(starts_with("Lineage"), names_to = "Curve", values_to = "weights")
#   )
# 
# df_w <- df_new %>%
#   dplyr::select(-pst) %>%
#   group_by(cells) %>%
#   mutate(weights = weights / sum(weights)) %>%
#   pivot_wider(names_from = "Curve", values_from = "weights")
# p6 <- ggplot(df_w, aes(x = Lineage1, y = Lineage3)) +
#   geom_hex() +
#   scale_fill_viridis_c(direction = -1) +
#   facet_wrap(~condition, scales = "free") +
#   geom_abline(slope = -1, intercept = 1, linetype = "dotted") +
#   geom_abline(slope = -1, intercept = 2/3, linetype = "dotted") +
#   geom_abline(slope = -1, intercept = 1/3, linetype = "dotted") +
#   annotate("text", x = .53, y = .53, label = "w3 = 0", angle = -52) +
#   annotate("text", x = .62, y = .1, label = "w3 = 1/3", angle = -52) +
#   annotate("text", x = .14, y = .14, label = "w3 = 2/3", angle = -52) +
#   theme(legend.position = "bottom") +
#   labs(x = "Weights for Lineage 1 (w1)", y = "Weights for Lineage 3 (w3)",
#        fill = "counts per hexagon")
# p6
# 
# df_w <- df_new %>%
#   dplyr::select(-pst) %>%
#   group_by(cells) %>%
#   mutate(weights = weights / sum(weights)) %>%
#   ungroup() %>%
#   group_by(condition, Curve) %>%
#   summarise(weights = mean(weights), .groups = NULL)
# 
# p7 <- ggplot(df_w, aes(x = Curve, fill = condition, y = weights)) +
#   geom_col(position = "dodge") +
#   scale_fill_brewer(palette = "Accent") +
#   theme(legend.position = c(.7, .7)) +
#   labs(x = "", y = "Mean weight")
# p7
# 
# write.csv(slingPseudotime(sds), file = '../output/monocle_2/cd4_memory/tradeSeq_all.csv')
# save(data, sds, sdss, file = '../output/monocle_2/cd4_memory/tradeSeq_all.RData')
load('../output/monocle_2/cd4_memory/tradeSeq_all.RData')
weights <- condiments:::.sling_reassign(sds)

set.seed(3)
rna <- seurat@assays$RNA@counts
filter <- apply(rna, 1, function(g) {
  sum(g >= 5) >= 15#10
})
rna_sub <- rna[filter, ]
# rna_sub <- rna_sub[,cellnames(sds)]

library(BiocParallel)
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 64
# icMat <- evaluateK(counts = rna_sub,
#                    pseudotime = slingPseudotime(sds, na = FALSE),
#                    cellWeights = weights,
#                    conditions = factor(data[,cellnames(sds)]$condition),
#                    nGenes = 100,
#                    k = 3:7,
#                    parallel = TRUE,
#                    BPPARAM = BPPARAM)

set.seed(3)
table(colnames(data) == rownames(slingPseudotime(sds, na = FALSE)))

# data@int_metadata$slingshot <- sds
table(colnames(data) == colnames(rna_sub))
data_new <- fitGAM(counts = rna_sub[, rownames(slingPseudotime(sds, na = FALSE))],
                   pseudotime = slingPseudotime(sds, na = FALSE),
                   cellWeights = weights,
                   conditions = colData(data)[rownames(slingPseudotime(sds, na = FALSE)),]$condition,
                   parallel = TRUE,
                   BPPARAM = BPPARAM,
                   nknots = 6)

save(data_new, file = '../output/monocle_2/cd4_memory/tradeSeq_data_new_all.RData')

condRes <- conditionTest(data_new, l2fc = log2(2), lineages = TRUE)

save(condRes, file = '../output/monocle_2/cd4_memory/tradeSeq_condRes_all.RData')

condRes$padj <- p.adjust(condRes$pvalue, "fdr")
condRes$padj_lineage1 <- p.adjust(condRes$pvalue_lineage1, "fdr")
condRes$padj_lineage2 <- p.adjust(condRes$pvalue_lineage2, "fdr")
condRes$padj_lineage3 <- p.adjust(condRes$pvalue_lineage3, "fdr")
condRes$padj_lineage4 <- p.adjust(condRes$pvalue_lineage4, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)

sum(condRes$padj <= 0.05, na.rm = TRUE)
sum(condRes$padj_lineage1 <= 0.05, na.rm = TRUE)
sum(condRes$padj_lineage2 <= 0.05, na.rm = TRUE)
sum(condRes$padj_lineage3 <= 0.05, na.rm = TRUE)
sum(condRes$padj_lineage4 <= 0.05, na.rm = TRUE)

conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]

library(ggsci)
scales <- c(pal_npg("nrc")(8))

oo <- order(condRes$waldStat, decreasing = TRUE)

# most significant gene
p8 <- plotSmoothers(data_new, counts(data_new),
                    gene = rownames(assays(data_new)$counts)[oo[1]],
                    alpha = 1, curvesCols = scales, sample = .3) +
  scale_color_manual(values = scales) +
  ggtitle(rownames(assays(data_new)$counts)[oo[1]])

p9 <- plotSmoothers(data_new, counts(data_new),
                    gene = 'IL17A',
                    alpha = 1, curvesCols = scales, sample = .3) +
  scale_color_manual(values = scales) +
  ggtitle('IL17A')

data_plot <- p9$data
data_plot <- data_plot[data_plot$lineage %in% c('Lineage 1_Act_IL1B_IL23', 'Lineage 1_Act_IL1B_IL23_PGE2'),]
ggplot(data_plot, aes(x = time, y = log(1 + gene_count), col = lineage))+
  geom_point()

### based on mean smoother
condRes$padj_lineage1 <- p.adjust(condRes$pvalue_lineage1, "fdr")
conditionGenes_lineage1 <- rownames(condRes)[condRes$padj_lineage1 <= 0.05]
conditionGenes_lineage1 <- conditionGenes_lineage1[!is.na(conditionGenes_lineage1)]

plot_heatmap <- function(gene, num, name){
  yhatSmooth <- predictSmooth(data_new, gene = gene, nPoints = 50, tidy = FALSE) %>%
    log1p()
  yhatSmoothScaled <- t(apply(yhatSmooth[, c(1:50+50*(num-1), 201:250+50*(num-1))], 1, scales::rescale))
  heatSmooth_Act_IL1B_IL23 <- pheatmap(yhatSmoothScaled[, 1:50],
                                       cluster_cols = FALSE,
                                       show_rownames = FALSE, show_colnames = FALSE, main = "Act_IL1B_IL23", legend = FALSE,
                                       silent = TRUE
  )
  
  matchingHeatmap_Act_IL1B_IL23_PGE2 <- pheatmap(yhatSmoothScaled[heatSmooth_Act_IL1B_IL23$tree_row$order, 51:100],
                                                 cluster_cols = FALSE, cluster_rows = FALSE,
                                                 show_rownames = FALSE, show_colnames = FALSE, main = "Act_IL1B_IL23_PGE2",
                                                 legend = FALSE, silent = TRUE
  )
  
  p <- plot_grid(heatSmooth_Act_IL1B_IL23[[4]], matchingHeatmap_Act_IL1B_IL23_PGE2[[4]], 
                 NULL, NULL, ncol = 2, rel_widths = c(1.4, 1.2), rel_heights = c(10, 1)) +
    draw_text(name, x = .5, y = .05)
  return(p)
}

plot_heatmap(conditionGenes_lineage1, 1, 'Lineage 1')
plot_heatmap(c('CREM', 'PTGER4', 'IL23R', 'RORA', 'IL17A', 'IL17F', 'TBX21', 'IFNG', 'LTA'), 1, 'Lineage 1')
