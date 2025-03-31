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
library(SummarizedExperiment)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

data <- readRDS('../output/tcell_annotated_updated_2_conditions.RDS')
data$celltype <- data$celltype_updated

data@reductions$umap <- data@reductions$atac.umap
data@reductions$umap@key <- 'UMAP_'

cds <- as.cell_data_set(data, assay = 'peaks')
rm(data)

cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

get_earliest_principal_node <- function(cds, celltype=celltype){
  cell_ids <- which(colData(cds)[, "celltype"] == celltype)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds, 'CD4+ Naive (Resting)'))

plot_cells(cds, color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE)

saveRDS(cds, file = '../output/monocle_2/cds_original.RDS')

p1 <- plot_cells(cds, color_cells_by = "celltype",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)
p2 <- plot_cells(cds, color_cells_by = "Activated1",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)
p3 <- plot_cells(cds, color_cells_by = "Resting1",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)
p4 <- plot_cells(cds, color_cells_by = "pseudotime",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)
p5 <- plot_cells(cds, color_cells_by = "sample",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)
p6 <- plot_cells(cds, color_cells_by = "condition",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)

pdf('../plots/monocle_2/UMAP_original.pdf', width = 24, height = 16)
wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3)
dev.off()

library(cowplot)
pdf('../plots/monocle_2/UMAP_original_lineage.pdf', width = 4, height = 4)
p1+theme_cowplot()+theme(legend.position = 'none')
dev.off()

#######
data <- readRDS('../output/tcell_annotated_updated_2_conditions.RDS')
data$celltype <- data$celltype_updated
cds <- new_cell_data_set(data@assays$peaks@counts,
                         cell_metadata = data@meta.data)

rm(data)
## Step 1: Normalize and pre-process the data
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 100, cores = 16, method = 'LSI')

## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "sample", cores = 16, preprocess_method = 'LSI')

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds, cores = 16)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)

## Step 5: Learn a graph
cds <- learn_graph(cds)

## Step 6: Order cells
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, celltype=celltype){
  cell_ids <- which(colData(cds)[, "celltype"] == celltype)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds, 'CD4+ Naive (Resting)'))

plot_cells(cds, color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE)

saveRDS(cds, file = '../output/monocle_2/cds.RDS')

p1 <- plot_cells(cds, color_cells_by = "celltype",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)
p2 <- plot_cells(cds, color_cells_by = "Activated1",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)
p3 <- plot_cells(cds, color_cells_by = "Resting1",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)
p4 <- plot_cells(cds, color_cells_by = "pseudotime",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)
p5 <- plot_cells(cds, color_cells_by = "sample",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)
p6 <- plot_cells(cds, color_cells_by = "condition",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE)

pdf('../plots/monocle_2/UMAP.pdf', width = 24, height = 16)
wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3)
dev.off()