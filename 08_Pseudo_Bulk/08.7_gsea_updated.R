library(clusterProfiler)
library(enrichplot)
library(data.table)
library(parallel)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(stringr)
library(viridis)
library(org.Hs.eg.db)
library(parallel)
library(pheatmap)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

rna <- fread('../output/pseudo_bulk_updated/DE_rna_DESeq2_paired.txt')

levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
           'CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
           'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1',
           'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17',
           'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh',
           'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other',
           'CD8+ Naive (Resting)', 'CD8+ Naive (Activated)',
           'CD8+ Regulatory',
           'CD8+ Memory (Resting)', 'CD8+ Memory (Activated)',
           'MAITs (Resting)', 'MAITs (Activated)', 'Gamma Delta'
)

result_combine_ps <- function(data){
  results = mutate(data, Expression = ifelse(data$padj < 0.05, "DEGs", "Not Significant"))
  results$Expression[is.na(results$Expression)] <- "Not Significant"
  results$Expression[which((results$log2FoldChange > 0) & (results$padj < 0.05))] <- "Up-regulated"
  results$Expression[which((results$log2FoldChange < 0) & (results$padj < 0.05))] <- "Down-regulated"
  # results$Expression[!results$gene %in% data2[data2$p_val_adj < 0.05,]$gene] <- "Not Significant"
  # p = ggplot(results, aes(avg_log2FC, -log10(p_val_adj))) + geom_point(aes(col = Expression)) +
  #   scale_color_manual(values = c("dodgerblue3", "darkgrey","firebrick")) +
  #   theme_bw(base_size = 12) +
  #   theme(text = element_text(size=16)) +
  #   theme(legend.position = "bottom") + geom_hline(yintercept = -log10(0.05), color = "grey1", lty = 2, lwd = 0.5)+
  #   ylab('-log10(FDR-adjusted P-value)')+
  #   xlab('log2(Fold Change)')+
  #   labs(col = '')+
  #   ggtitle(title)+
  #   theme(plot.title = element_text(hjust = 0.5))
  return(results)
}

rna_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  print(celltype)
  data <- rna[rna$cell_type == celltype,]
  rna_list[[celltype]] <- result_combine_ps(data)
}

bulk <- fread('~/RWorkSpace/Duerr_bulk/RNA/gene/analysis_new/results_Act_IL1B_IL23_PGE2_vs_Act_IL1B_IL23.txt')
bulk <- tibble::column_to_rownames(bulk, var = 'V1')
bulk = mutate(bulk, Expression = ifelse(bulk$padj < 0.05, "DEGs", "Not Significant"))
bulk$Expression[is.na(bulk$Expression)] <- "Not Significant"
bulk$Expression[which((bulk$log2FoldChange > 0) & (bulk$padj < 0.05))] <- "Up-regulated"
bulk$Expression[which((bulk$log2FoldChange < 0) & (bulk$padj < 0.05))] <- "Down-regulated"

draw_plot <- function(df){
  # we want the log2 fold change 
  original_gene_list <- df$log2FoldChange
  
  # name the vector
  names(original_gene_list) <- df$gene
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  gse <- gseGO(geneList=gene_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "BH")
  
  gse_m <- gse@result
  gse_m <- gse_m[order(gse_m$NES),]
  # gse_m_filt = rbind(head(gse_m, n = 5),
  #                    tail(gse_m, n = 5))
  # gse_m_filt <- gse_m_filt[order(gse_m_filt$NES, decreasing = F),]
  # gse_m_filt$Description <- factor(gse_m_filt$Description, levels = gse_m_filt$Description)
  # gse_m_filt$Description_re <- factor(str_wrap(str_to_sentence(gse_m_filt$Description), 40), levels = str_wrap(str_to_sentence(gse_m_filt$Description), 40))
  # gse_m_filt$log10padj <- -log10(gse_m_filt$p.adjust)
  # p <- ggplot(gse_m_filt, aes(Description_re, NES)) +
  #   geom_segment(aes(xend=Description_re, y=0, yend=NES)) +
  #   geom_point(aes( fill = log10padj, size = setSize),
  #              shape=21) +
  #   scale_fill_viridis(direction = -1, option = 'C') +
  #   coord_flip() +
  #   labs(x="", y="Normalized Enrichment Score", fill = '-log10(FDR)',
  #        title="GSEA - Biological Processes", subtitle = subtitle) + 
  #   theme_cowplot()+
  #   theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(hjust = 0))
  return(gse_m)
}

p_list = mclapply(levels[1:20], function(i){
  p = draw_plot(rna_list[[i]])
  return(p)
}, mc.cores = 4)
names(p_list) = levels

p_list[['Bulk']] = draw_plot(bulk)

p_list_re = lapply(1:21, function(i){
  celltype = names(p_list)[i]
  df = cbind(celltype, p_list[[i]])
  df$p_val_adj = p.adjust(df$pvalue, method = 'fdr')
  df$log10_p_val_adj = -log10(df$p_val_adj)
  return(df)
})
names(p_list_re) = names(p_list)

gsea = Reduce(rbind, p_list_re)
gsea$celltype = factor(gsea$celltype, levels = c(levels, 'Bulk'))

mat = dcast(gsea[,c(1,3,14)], Description ~ celltype)
mat = tibble::column_to_rownames(mat, 'Description')
mat = as.matrix(mat)
mat[is.na(mat)] = 1
table(rowMax(mat) > -log10(0.00000001))
mat = mat[rowMax(mat) > -log10(0.00000001),]

p <- pheatmap(mat, cluster_rows = T, cluster_cols = F)
p

order <- p$tree_row$labels[p$tree_row$order]

mat1 = dcast(gsea[,c(1,3,6)], Description ~ celltype)
mat1 = tibble::column_to_rownames(mat1, 'Description')
mat1 = as.matrix(mat1)
mat1[is.na(mat1)] = 0
mat1 = mat1[order,]

p1 <- pheatmap(mat1, cluster_rows = T, cluster_cols = F)
p1

order1 <- p1$tree_row$labels[p1$tree_row$order]

index = which(gsea$Description %in% order)
gsea_sub = gsea[index,]
gsea_sub$Description <- factor(gsea_sub$Description, levels = order)

pdf("../plots/pseudo_bulk/gsea/gsea_fdr_0.00000001.pdf", width = 12, height = 24)
ggplot(gsea_sub, aes(x = celltype, y = Description, size = log10_p_val_adj, fill = NES))+
    geom_point(shape=21)+
    scale_fill_gradient2(low = 'dodgerblue3', mid = 'white', high = 'firebrick', na.value = 'darkgray')+
    labs(x = '', y = '', size = '-log10(FDR)')+
    theme_cowplot()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

gsea_sub$Description <- factor(gsea_sub$Description, levels = order1)

pdf("../plots/pseudo_bulk/gsea/gsea_fdr_0.00000001_re.pdf", width = 12, height = 24)
ggplot(gsea_sub, aes(x = celltype, y = Description, size = log10_p_val_adj, fill = NES))+
    geom_point(shape=21)+
    scale_fill_gradient2(low = 'dodgerblue3', mid = 'white', high = 'firebrick', na.value = 'darkgray')+
    labs(x = '', y = '', size = '-log10(FDR)')+
    theme_cowplot()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()
