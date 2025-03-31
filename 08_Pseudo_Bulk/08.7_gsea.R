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

draw_plot <- function(df, subtitle){
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
  gse_m_filt = rbind(head(gse_m, n = 5),
                     tail(gse_m, n = 5))
  gse_m_filt <- gse_m_filt[order(gse_m_filt$NES, decreasing = F),]
  gse_m_filt$Description <- factor(gse_m_filt$Description, levels = gse_m_filt$Description)
  gse_m_filt$Description_re <- factor(str_wrap(str_to_sentence(gse_m_filt$Description), 40), levels = str_wrap(str_to_sentence(gse_m_filt$Description), 40))
  gse_m_filt$log10padj <- -log10(gse_m_filt$p.adjust)
  p <- ggplot(gse_m_filt, aes(Description_re, NES)) +
    geom_segment(aes(xend=Description_re, y=0, yend=NES)) +
    geom_point(aes( fill = log10padj, size = setSize),
               shape=21) +
    scale_fill_viridis(direction = -1, option = 'C') +
    coord_flip() +
    labs(x="", y="Normalized Enrichment Score", fill = '-log10(FDR)',
         title="GSEA - Biological Processes", subtitle = subtitle) + 
    theme_cowplot()+
    theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(hjust = 0))
  return(p)
}

p_list = mclapply(1:20, function(i){
  p = draw_plot(rna_list[[levels[i]]], levels[i])
  return(p)
}, mc.cores = 4)

for (i in 1:20){
  assign(paste0('p', i), p_list[[i]])
}

p21 = draw_plot(bulk, 'Bulk')

pdf('../plots/pseudo_bulk/gsea/gsea.pdf', width = 21, height = 42)
wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21, ncol = 3)
dev.off()
