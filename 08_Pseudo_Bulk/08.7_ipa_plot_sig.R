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
library(pheatmap)

setwd("~/RWorkSpace/DOGMA-seq/PBMC/code")

levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
           'CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
           'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1',
           'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17',
           'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh',
           'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other',
           'CD8+ Naive (Resting)', 'CD8+ Naive (Activated)',
           'CD8+ Regulatory',
           'CD8+ Memory (Resting)', 'CD8+ Memory (Activated)',
           'MAITs (Resting)', 'MAITs (Activated)', 'Gamma Delta', 'Bulk'
)

ipa_list = lapply(levels[1:21], function(celltype){
    df = fread(paste0("../plots/pseudo_bulk/RNA_updated/IPA_results/", celltype, ".txt"))[,-6]
    print(nrow(df))
    df$p_val = exp(-df$`-log(p-value)`)
    df$p_val_adj = p.adjust(df$p_val, method = 'fdr')
    df$log10_p_val_adj = -log10(df$p_val_adj)
    df = cbind(celltype, df)
    return(df)
})
ipa = Reduce(rbind, ipa_list)
ipa$celltype = factor(ipa$celltype, levels = levels)

mat = dcast(ipa[,c(1:2,9)], `Ingenuity Canonical Pathways` ~ celltype)
mat = tibble::column_to_rownames(mat, 'Ingenuity Canonical Pathways')
mat = as.matrix(mat)
mat[is.na(mat)] = 1
table(rowMax(mat) > -log10(0.0001))
mat = mat[rowMax(mat) > -log10(0.0001),]

p <- pheatmap(mat, cluster_rows = T, cluster_cols = F)
p

order <- p$tree_row$labels[p$tree_row$order]

mat1 = dcast(ipa[,c(1:2,5)], `Ingenuity Canonical Pathways` ~ celltype)
mat1 = tibble::column_to_rownames(mat1, 'Ingenuity Canonical Pathways')
mat1 = as.matrix(mat1)
mat1[is.na(mat1)] = 0
mat1 = mat1[order,]

p1 <- pheatmap(mat1, cluster_rows = T, cluster_cols = F)
p1

order1 <- p1$tree_row$labels[p1$tree_row$order]

index = which(ipa$`Ingenuity Canonical Pathways` %in% order)
ipa_sub = ipa[index,]
ipa_sub$log10_p_val_adj[!ipa_sub$p_val < 0.05] = NA
ipa_sub$`Ingenuity Canonical Pathways` <- factor(ipa_sub$`Ingenuity Canonical Pathways`, levels = order)

pdf("../plots/pseudo_bulk/ipa/ipa_fdr_0.0001_sig.pdf", width = 10, height = 10)
ggplot(ipa_sub, aes(x = celltype, y = `Ingenuity Canonical Pathways`, size = log10_p_val_adj, fill = `z-score`))+
    geom_point(shape=21)+
    scale_fill_gradient2(low = 'dodgerblue3', mid = 'white', high = 'firebrick', na.value = 'darkgray')+
    labs(x = '', y = '', size = '-log10(FDR)')+
    theme_cowplot()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

ipa_sub$`Ingenuity Canonical Pathways` <- factor(ipa_sub$`Ingenuity Canonical Pathways`, levels = order1)

pdf("../plots/pseudo_bulk/ipa/ipa_fdr_0.0001_re_sig.pdf", width = 10, height = 10)
ggplot(ipa_sub, aes(x = celltype, y = `Ingenuity Canonical Pathways`, size = log10_p_val_adj, fill = `z-score`))+
    geom_point(shape=21)+
    scale_fill_gradient2(low = 'dodgerblue3', mid = 'white', high = 'firebrick', na.value = 'darkgray')+
    labs(x = '', y = '', size = '-log10(FDR)')+
    theme_cowplot()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

order2 = order1[c(1:9,19:25,10:18,26:32,38:39,33:37)]

ipa_sub$`Ingenuity Canonical Pathways` <- factor(ipa_sub$`Ingenuity Canonical Pathways`, levels = order2)

pdf("../plots/pseudo_bulk/ipa/ipa_fdr_0.0001_re_re_sig.pdf", width = 10, height = 10)
ggplot(ipa_sub, aes(x = celltype, y = `Ingenuity Canonical Pathways`, size = log10_p_val_adj, fill = `z-score`))+
    geom_point(shape=21)+
    scale_fill_gradient2(low = 'dodgerblue3', mid = 'white', high = 'firebrick', na.value = 'darkgray')+
    labs(x = '', y = '', size = '-log10(FDR)')+
    theme_cowplot()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

View(ipa_sub[str_detect(ipa_sub$Molecules, 'CREM'),])
table(as.character(ipa_sub[str_detect(ipa_sub$Molecules, 'CREM'),]$`Ingenuity Canonical Pathways`))
# Protein Kinase A Signaling 
#                         20 
View(ipa_sub[str_detect(ipa_sub$Molecules, 'STAT4'),])
table(as.character(ipa_sub[str_detect(ipa_sub$Molecules, 'STAT4'),]$`Ingenuity Canonical Pathways`))
# IL-12 Signaling and Production in Macrophages 
#                                            20 
#                       IL-23 Signaling Pathway 
#                                            20 
#                            JAK/STAT Signaling 
#                                            20 
#                 Natural Killer Cell Signaling 
#                                            20 
#                Th1 and Th2 Activation Pathway 
#                                            20
View(ipa_sub[str_detect(ipa_sub$Molecules, 'PDE4D'),])
table(as.character(ipa_sub[str_detect(ipa_sub$Molecules, 'PDE4D'),]$`Ingenuity Canonical Pathways`))
# Protein Kinase A Signaling 
#                         20 
View(ipa[str_detect(ipa$`Ingenuity Canonical Pathways`, 'EIF2'),])
################
mat = dcast(ipa[,c(1:2,9)], `Ingenuity Canonical Pathways` ~ celltype)
mat = tibble::column_to_rownames(mat, 'Ingenuity Canonical Pathways')
mat = as.matrix(mat)
mat[is.na(mat)] = 1
table(rowMax(mat) > -log10(0.001))
mat = mat[rowMax(mat) > -log10(0.001),]

p <- pheatmap(mat, cluster_rows = T, cluster_cols = F)
p

order <- p$tree_row$labels[p$tree_row$order]

mat1 = dcast(ipa[,c(1:2,5)], `Ingenuity Canonical Pathways` ~ celltype)
mat1 = tibble::column_to_rownames(mat1, 'Ingenuity Canonical Pathways')
mat1 = as.matrix(mat1)
mat1[is.na(mat1)] = 0
mat1 = mat1[order,]

p1 <- pheatmap(mat1, cluster_rows = T, cluster_cols = F)
p1

order1 <- p1$tree_row$labels[p1$tree_row$order]

index = which(ipa$`Ingenuity Canonical Pathways` %in% order)
ipa_sub = ipa[index,]
ipa_sub$log10_p_val_adj[!ipa_sub$p_val < 0.05] = NA
ipa_sub$`Ingenuity Canonical Pathways` <- factor(ipa_sub$`Ingenuity Canonical Pathways`, levels = order)

pdf("../plots/pseudo_bulk/ipa/ipa_fdr_0.001_sig.pdf", width = 12, height = 18)
ggplot(ipa_sub, aes(x = celltype, y = `Ingenuity Canonical Pathways`, size = log10_p_val_adj, fill = `z-score`))+
    geom_point(shape=21)+
    scale_fill_gradient2(low = 'dodgerblue3', mid = 'white', high = 'firebrick', na.value = 'darkgray')+
    labs(x = '', y = '', size = '-log10(FDR)')+
    theme_cowplot()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

ipa_sub$`Ingenuity Canonical Pathways` <- factor(ipa_sub$`Ingenuity Canonical Pathways`, levels = order1)

pdf("../plots/pseudo_bulk/ipa/ipa_fdr_0.001_re_sig.pdf", width = 12, height = 18)
ggplot(ipa_sub, aes(x = celltype, y = `Ingenuity Canonical Pathways`, size = log10_p_val_adj, fill = `z-score`))+
    geom_point(shape=21)+
    scale_fill_gradient2(low = 'dodgerblue3', mid = 'white', high = 'firebrick', na.value = 'darkgray')+
    labs(x = '', y = '', size = '-log10(FDR)')+
    theme_cowplot()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()
