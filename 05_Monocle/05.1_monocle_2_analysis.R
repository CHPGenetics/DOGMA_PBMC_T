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
library(ComplexHeatmap)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

# cds <- readRDS('../output/monocle_2/cds_original.RDS')
# 
# get_earliest_principal_node <- function(cds, celltype=celltype){
#   cell_ids <- which(colData(cds)[, "celltype"] == celltype)
#   
#   closest_vertex <-
#     cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#   root_pr_nodes <-
#     igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
#                                                               (which.max(table(closest_vertex[cell_ids,]))))]
#   
#   root_pr_nodes
# }
# 
# my_get_earliest_principal_node <- function(cds, celltype=celltype){
#   cell_ids <- which(colData(cds)[, "celltype"] %in% celltype)
#   
#   closest_vertex <-
#     cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#   # root_pr_nodes <-
#   #   igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
#   #                                                             (which.max(table(closest_vertex[cell_ids,]))))]
#   #
#   # root_pr_nodes
#   df <- as.data.frame(table(closest_vertex[cell_ids,]))
#   return(df)
# }
# 
# # data <- readRDS('../output/all_qc_updated/tcell.RDS')
# #
# # chromvar <- data@assays$chromvar@data
# cds_sub <- cds[,which(colData(cds)[, "celltype"] %in% c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)'))]
# df <- my_get_earliest_principal_node(cds_sub, c('CD4+ Naive (Resting)'))
# cds_sub <- order_cells(cds_sub, root_pr_nodes='Y_539')
# 
# p1 <- plot_cells(cds_sub, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# # cells <- read.csv('../output/monocle/cd4_naive.csv', row.names = 'X')
# # table(colnames(cds_sub) %in% rownames(cells))
# # cds_sub_cut <- cds_sub[,which(colnames(cds_sub) %in% rownames(cells))]
# pseudotime(cds_sub)
# cds_sub_cut <- choose_cells(cds_sub)
# max(pseudotime(cds_sub_cut))
# cds_sub_cut <- cds_sub[,which(pseudotime(cds_sub) < 14.35621)]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# p2 <- plot_cells(cds_sub_cut, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# 
# p1 | p2
# df <- my_get_earliest_principal_node(cds_sub_cut, 'CD4+ Naive (Resting)')
# cds_sub_cut <- order_cells(cds_sub_cut, root_pr_nodes='Y_539')
# p3 <- plot_cells(cds_sub_cut, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# p1 | p2 | p3
# save(cds_sub, cds_sub_cut, file = '../output/monocle_2/cd4_naive.RData')
# 
# write.csv(pseudotime(cds_sub_cut), file = '../output/monocle_2/cd4_naive.csv')
# 
# ########
# cds_sub <- cds[,which(colData(cds)[, "celltype"] %in% c('CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1',
#                                                         'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17',
#                                                         'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh',
#                                                         'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other'))]
# df <- my_get_earliest_principal_node(cds_sub, c('CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Resting) - Other'))
# cds_sub <- order_cells(cds_sub, root_pr_nodes='Y_561')
# 
# p1 <- plot_cells(cds_sub, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# # cells <- read.csv('../output/monocle/cd4_memory.csv', row.names = 'X')
# # table(colnames(cds_sub) %in% rownames(cells))
# # cds_sub_cut <- cds_sub[,which(colnames(cds_sub) %in% rownames(cells))]
# pseudotime(cds_sub)
# cds_sub_cut <- choose_cells(cds_sub)
# max(pseudotime(cds_sub_cut))
# cds_sub_cut <- cds_sub[,which(pseudotime(cds_sub) < 16)]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# p2 <- plot_cells(cds_sub_cut, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# 
# p1 | p2
# save(cds_sub, cds_sub_cut, file = '../output/monocle_2/cd4_memory.RData')
# 
# write.csv(pseudotime(cds_sub_cut), file = '../output/monocle_2/cd4_memory.csv')
# 
# ########
# cds_sub <- cds[,which(colData(cds)[, "celltype"] %in% c('CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)'))]
# df <- my_get_earliest_principal_node(cds_sub, c('CD4+ Regulatory (Resting)'))
# cds_sub <- order_cells(cds_sub, root_pr_nodes='Y_380')
# 
# p1 <- plot_cells(cds_sub, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# # cells <- read.csv('../output/monocle/cd4_memory.csv', row.names = 'X')
# # table(colnames(cds_sub) %in% rownames(cells))
# # cds_sub_cut <- cds_sub[,which(colnames(cds_sub) %in% rownames(cells))]
# pseudotime(cds_sub)
# cds_sub_cut <- choose_cells(cds_sub)
# max(pseudotime(cds_sub_cut))
# cds_sub_cut <- cds_sub[,which(pseudotime(cds_sub) < 4.69144)]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# p2 <- plot_cells(cds_sub_cut, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# 
# p1 | p2
# save(cds_sub, cds_sub_cut, file = '../output/monocle_2/cd4_regulatory.RData')
# 
# write.csv(pseudotime(cds_sub_cut), file = '../output/monocle_2/cd4_regulatory.csv')
# 
# ########
# cds_sub <- cds[,which(colData(cds)[, "celltype"] %in% c('CD8+ Naive (Resting)', 'CD8+ Naive (Activated)'))]
# df <- my_get_earliest_principal_node(cds_sub, c('CD8+ Naive (Resting)'))
# cds_sub <- order_cells(cds_sub, root_pr_nodes='Y_384')
# 
# p1 <- plot_cells(cds_sub, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# # cells <- read.csv('../output/monocle/cd8_naive.csv', row.names = 'X')
# # table(colnames(cds_sub) %in% rownames(cells))
# # cds_sub_cut <- cds_sub[,which(colnames(cds_sub) %in% rownames(cells))]
# pseudotime(cds_sub)
# cds_sub_cut <- choose_cells(cds_sub)
# max(pseudotime(cds_sub_cut))
# cds_sub_cut <- cds_sub[,which(pseudotime(cds_sub) < 9.45313)]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_cells(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# p2 <- plot_cells(cds_sub_cut, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# 
# p1 | p2
# df <- my_get_earliest_principal_node(cds_sub_cut, 'CD8+ Naive (Resting)')
# cds_sub_cut <- order_cells(cds_sub_cut, root_pr_nodes='Y_733')
# p3 <- plot_cells(cds_sub_cut, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# p1 | p2 | p3
# save(cds_sub, cds_sub_cut, file = '../output/monocle_2/cd8_naive.RData')
# 
# write.csv(pseudotime(cds_sub_cut), file = '../output/monocle_2/cd8_naive.csv')
# 
# ########
# cds_sub <- cds[,which(colData(cds)[, "celltype"] %in% c('CD8+ Memory (Resting)', 'CD8+ Memory (Activated)'))]
# df <- my_get_earliest_principal_node(cds_sub, c('CD8+ Memory (Resting)'))
# cds_sub <- order_cells(cds_sub, root_pr_nodes='Y_898')
# 
# p1 <- plot_cells(cds_sub, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# pseudotime(cds_sub)
# cds_sub_cut <- choose_cells(cds_sub)
# max(pseudotime(cds_sub_cut))
# cds_sub_cut <- cds_sub[,which(pseudotime(cds_sub) < 13.57058)]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# p2 <- plot_cells(cds_sub_cut, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# 
# p1 | p2
# save(cds_sub, cds_sub_cut, file = '../output/monocle_2/cd8_memory.RData')
# 
# write.csv(pseudotime(cds_sub_cut), file = '../output/monocle_2/cd8_memory.csv')
# 
# ########
# cds_sub <- cds[,which(colData(cds)[, "celltype"] %in% c('MAITs (Resting)', 'MAITs (Activated)'))]
# 
# df <- my_get_earliest_principal_node(cds_sub, c('MAITs (Resting)'))
# cds_sub <- order_cells(cds_sub, root_pr_nodes='Y_394')
# 
# p1 <- plot_cells(cds_sub, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# pseudotime(cds_sub)
# cds_sub_cut <- choose_cells(cds_sub)
# max(pseudotime(cds_sub_cut))
# cds_sub_cut <- cds_sub[,which(pseudotime(cds_sub) < 6.070789)]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# cds_sub_cut_cut <- choose_graph_segments(cds_sub_cut)
# cds_sub_cut <- cds_sub_cut[,which(!colnames(cds_sub_cut) %in% colnames(cds_sub_cut_cut))]
# p2 <- plot_cells(cds_sub_cut, color_cells_by = "pseudotime",
#                  label_groups_by_cluster=FALSE,
#                  label_leaves=FALSE,
#                  label_branch_points=FALSE)
# 
# p1 | p2
# save(cds_sub, cds_sub_cut, file = '../output/monocle_2/mait.RData')
# 
# write.csv(pseudotime(cds_sub_cut), file = '../output/monocle_2/mait.csv')

#####################
data <- readRDS('../output/tcell_annotated_updated_2_conditions.RDS')
add_pseudotime <- function(name){
  pseudotime <- read.csv(paste0('../output/monocle_2/', name, '.csv'), row.names = 'X')
  data <- AddMetaData(data, pseudotime, paste0('monocle_', name))
  return(data)
}
data <- add_pseudotime('cd4_naive')
data <- add_pseudotime('cd4_memory')
data <- add_pseudotime('cd4_regulatory')
data <- add_pseudotime('cd8_naive')
data <- add_pseudotime('cd8_memory')
data <- add_pseudotime('mait')

# pseudotime <- read.csv(paste0('../output/monocle/', 'cd4_memory/Th17_Th1_3', '.csv'), row.names = 'X')
# data <- AddMetaData(data, pseudotime, paste0('monocle_', 'Th17_Th1'))

plot_pseudotime <- function(name, title){
  p <- FeaturePlot(data, paste0('monocle_', name), reduction = 'atac.umap') + scale_color_viridis_c(direction = -1) + ggtitle(title)
  pdf(paste0('../plots/monocle_2/monocle_', name, '.pdf'), width = 4, height = 4)
  print(p)
  dev.off()
  return(p)
}

p1 <- plot_pseudotime('cd4_naive', 'CD4+ naive T cells')
p2 <- plot_pseudotime('cd4_memory', 'CD4+ memory T cells')
p3 <- plot_pseudotime('cd4_regulatory', 'CD4+ regulatory T cells')
p4 <- plot_pseudotime('cd8_naive', 'CD8+ naive T cells')
p5 <- plot_pseudotime('cd8_memory', 'CD8+ memory T cells')
p6 <- plot_pseudotime('mait', 'MAITs')
# plot_pseudotime('Th17_Th1', '')

pdf(paste0('../plots/monocle_2/monocle_plots.pdf'), width = 12, height = 8)
wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3)
dev.off()

# tradeseq <- read.csv('../output/monocle/cd4_memory/tradeSeq_all.csv', row.names = 'X')
# data <- AddMetaData(data, tradeseq)
# 
# p1 <- FeaturePlot(data, "Lineage1", reduction = 'atac.umap') + scale_color_viridis_c(direction = -1) + ggtitle('CD4+ memory T cells (Lineage 1)')
# p2 <- FeaturePlot(data, "Lineage2", reduction = 'atac.umap') + scale_color_viridis_c(direction = -1) + ggtitle('CD4+ memory T cells (Lineage 2)')
# p3 <- FeaturePlot(data, "Lineage3", reduction = 'atac.umap') + scale_color_viridis_c(direction = -1) + ggtitle('CD4+ memory T cells (Lineage 3)')
# p4 <- FeaturePlot(data, "Lineage4", reduction = 'atac.umap') + scale_color_viridis_c(direction = -1) + ggtitle('CD4+ memory T cells (Lineage 4)')
# pdf(paste0('../plots/monocle/tradeSeq_plots.pdf'), width = 16, height = 4)
# wrap_plots(p1, p2, p3, p4, ncol = 4)
# dev.off()
