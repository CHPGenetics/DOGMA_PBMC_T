library(data.table)
library(parallel)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(venn)
library(ggvenn)
library(ggsci)
library(ggrepel)
library(plyr)
library(doMC)
library(ArchR)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

ibd_peak <- read.table('../output/SNP/overlap_ibd_peak_snp.txt')

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

reform_p2g <- function(p2g){
  peak <- metadata(p2g)$peakSet
  gene <- metadata(p2g)$geneSet
  p2g_df <- as.data.frame(p2g)
  
  peak_sub <- peak[p2g_df$idxATAC]
  gene_sub <- gene[p2g_df$idxRNA]
  peak_df <- data.frame(chr = seqnames(peak_sub), start = start(peak_sub), end = end(peak_sub))
  gene_df <- data.frame(chr = seqnames(gene_sub), pos = start(gene_sub), gene = gene_sub$name)
  peak_df$chr_pos <- paste(peak_df$chr, peak_df$start, peak_df$end, sep = '-')
  
  p2g_df$peak_chr <- peak_df$chr
  p2g_df$peak_start <- peak_df$start
  p2g_df$peak_end <- peak_df$end
  p2g_df$peak_chr_pos <- peak_df$chr_pos
  p2g_df$gene_chr <- gene_df$chr
  p2g_df$gene_pos <- gene_df$pos
  p2g_df$gene <- gene_df$gene
  length(unique(p2g_df$peak_chr_pos))#188999
  length(unique(p2g_df$gene))#25736
  p2g_df$distance <- pmin(abs(p2g_df$peak_start - p2g_df$gene_pos), abs(p2g_df$peak_end - p2g_df$gene_pos))
  p2g_df$distance[p2g_df$peak_start <= p2g_df$gene_pos & p2g_df$gene_pos <= p2g_df$peak_end] <- 0
  # p2g_df$gene_merged <- NA
  # p2g_df$gene_pos_merged <- NA
  # p2g_df$gene_num <- NA
  # for (i in 1:nrow(p2g_df)){
  #   tmp <- p2g_df[p2g_df$peak_chr_pos == p2g_df$peak_chr_pos[i],]
  #   p2g_df$gene_merged[i] <- paste(tmp$gene, collapse = '|')
  #   p2g_df$gene_pos_merged[i] <- paste(tmp$gene_pos, collapse = '|')
  #   p2g_df$gene_num[i] <- nrow(tmp)
  # }
  
  # p2g_df_dedup <- p2g_df[!duplicated(p2g_df$peak_chr_pos),]
  # p2g_df_dedup <- p2g_df_dedup[,-c(11:14)]
  
  # merged <- ddply(p2g_df, 'peak_chr_pos', function(tmp){
  #   gene_merged <- paste(tmp$gene, collapse = '|')
  #   distance_merged <- paste(tmp$distance, collapse = '|')
  #   gene_num <- nrow(tmp)
  #   correlation_merged <- paste(tmp$Correlation, collapse = '|')
  #   fdr_merged <- paste(tmp$FDR, collapse = '|')
  #   return(c(gene_merged, distance_merged, gene_num, correlation_merged, fdr_merged))
  # }, .parallel = T)
  # merged <- mclapply(1:nrow(p2g_df_dedup), function(i){
  #   tmp <- p2g_df[p2g_df$peak_chr_pos == p2g_df_dedup$peak_chr_pos[i],]
  #   gene_merged <- paste(tmp$gene, collapse = '|')
  #   gene_pos_merged <- paste(tmp$gene_pos, collapse = '|')
  #   distance_merged <- paste(tmp$distance, collapse = '|')
  #   gene_num <- nrow(tmp)
  #   return(c(gene_merged, gene_pos_merged, distance_merged, gene_num))
  # }, mc.cores = 24)
  # 
  # merged_df <- as.data.frame(matrix(unlist(merged), ncol = 4, byrow = T))
  # colnames(merged_df) <- c('gene_merged', 'gene_pos_merged', 'distance_merged', 'gene_num')
  # colnames(merged)[2:6] <- c('gene_merged', 'distance_merged', 'gene_num', 'correlation_merged', 'fdr_merged')
  p2g_df <- p2g_df[,-c(1:2)]
  return(p2g_df)
}

p2g_celltype_link <- list()
for(i in 1:9){
  celltype <- levels[i]
  name <- paste0('0',i,'_',celltype)
  p2g <- readRDS(paste0('../output/ArchR/DOGMA_filtered_multiome/', name, '/p2g_fdr_0.001.RDS'))
  p2g_celltype_link[[celltype]] <- reform_p2g(p2g)
}
for(i in 10:20){
  celltype <- levels[i]
  name <- paste0(i,'_',celltype)
  p2g <- readRDS(paste0('../output/ArchR/DOGMA_filtered_multiome/', name, '/p2g_fdr_0.001.RDS'))
  p2g_celltype_link[[celltype]] <- reform_p2g(p2g)
}

check_pair = function(celltype){
  df <- p2g_celltype_link[[celltype]]
  table = data.frame(celltype = celltype,
    pos_pair = as.numeric(table(df$Correlation > 0)['TRUE']),
    neg_pair = as.numeric(table(df$Correlation > 0)['FALSE']),
    genes = length(unique(df$gene)),
    peaks = length(unique(df$peak_chr_pos)))
  return(table)
}

df_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  df <- check_pair(celltype)
  df_list[[celltype]] <- df
}
df <- Reduce(rbind, df_list)

data_plot = data.frame(celltype = df$celltype, pair = df$pos_pair, group = 'Positive')
data_plot = rbind(data_plot, data.frame(celltype = df$celltype, pair = df$neg_pair, group = 'Negative'))
data_plot$celltype = factor(data_plot$celltype, levels = levels)

library(scales)
p1 = ggplot(data_plot, aes(x = pair, y = celltype, fill = group))+
  geom_col(position = 'stack', col = 'white')+
  geom_text(aes(label = celltype, x = 2e4), hjust = 0, colour = "black")+
  theme_bw()+
  scale_fill_manual(values = c('Positive' = 'firebrick', 'Negative' = 'dodgerblue3'))+
  scale_x_continuous(labels = scientific, n.breaks = 13)+
  labs(y = 'Cell types', x = 'Peak-to-gene linkages within 1M bps (cis-eQTA)', fill = '')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'bottom', axis.text.y = element_blank(), axis.ticks.y = element_blank())
pdf('../plots/eQTA/pairs.pdf', width = 5, height = 5)
p1
dev.off()

df$celltype = factor(df$celltype, levels = levels)
p1 = ggplot(df, aes(x = genes, y = celltype, fill = celltype))+
  geom_col(position = 'stack', col = 'white')+
  geom_text(aes(label = celltype, x = 200), hjust = 0, colour = "black")+
  theme_bw()+
  scale_x_continuous(n.breaks = 14)+
  labs(y = 'Cell types', x = '# of eGenes', fill = '')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', axis.text.y = element_blank(), axis.ticks.y = element_blank())
pdf('../plots/eQTA/genes.pdf', width = 3, height = 5)
p1
dev.off()

p1 = ggplot(df, aes(x = peaks, y = celltype, fill = celltype))+
  geom_col(position = 'stack', col = 'white')+
  geom_text(aes(label = celltype, x = 2000), hjust = 0, colour = "black")+
  theme_bw()+
  scale_x_continuous(n.breaks = 10)+
  labs(y = 'Cell types', x = '# of ePeaks', fill = '')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', axis.text.y = element_blank(), axis.ticks.y = element_blank())
pdf('../plots/eQTA/peaks.pdf', width = 3, height = 5)
p1
dev.off()

###############
check_gene = function(celltype){
  df <- p2g_celltype_link[[celltype]]
  table = unique(df$gene)
  return(table)
}

df_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  df <- check_gene(celltype)
  df_list[[celltype]] <- df
}

check_overlap = function(celltype1, celltype2){
  ft1 = df_list[[celltype1]]
  ft2 = df_list[[celltype2]]
  ft_shared = intersect(ft1, ft2)
  return(length(ft_shared)/length(ft1))
}

mat = matrix(NA, 20, 20)
for (i in 1:20){
  for (j in 1:20){
    mat[i,j] = check_overlap(i,j)
  }
}
rownames(mat) = colnames(mat) = levels

library(ComplexHeatmap)
pdf('../plots/eQTA/shared_gene.pdf', width = 8, height = 5)
Heatmap(mat*100, cluster_rows = F, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = "%"), col = viridis(100, option = "C", direction = -1))
dev.off()

###############
check_peak = function(celltype){
  df <- p2g_celltype_link[[celltype]]
  table = unique(df$peak_chr_pos)
  return(table)
}

df_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  df <- check_peak(celltype)
  df_list[[celltype]] <- df
}

mat = matrix(NA, 20, 20)
for (i in 1:20){
  for (j in 1:20){
    mat[i,j] = check_overlap(i,j)
  }
}
rownames(mat) = colnames(mat) = levels

pdf('../plots/eQTA/shared_peak.pdf', width = 8, height = 5)
Heatmap(mat*100, cluster_rows = F, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = "%"), col = viridis(100, option = "C", direction = -1))
dev.off()

###############
#peak_annotation <- readRDS('../output/pseudo_bulk/peak_annotation.RDS')
#peak_annotation$name = rownames(peak_annotation)

anno <- fread('../output/pseudo_bulk/peaks.annotation_sorted.txt')[-1]
anno$name <- paste(anno$Chr, anno$Start-1, anno$End, sep = '-')

check_top_eQTA = function(celltype){
  tmp = p2g_celltype_link[[celltype]]
  tmp$name = tmp$peak_chr_pos

  tmp1=tmp[order(tmp$FDR),]
  tmp2=tmp1[!duplicated(tmp1$name),]
  tmp3=merge(tmp2, anno, by = 'name', all.x = T)
  table=data.frame(celltype = celltype,
  top_nearest = as.numeric(table(tmp3$gene == tmp3$`Gene Name`)['TRUE']),
  top_not_nearest = as.numeric(table(tmp3$gene == tmp3$`Gene Name`)['FALSE']))
  return(table)
}

top_eQTA_list = lapply(levels[1:20],function(i){check_top_eQTA(i)})
top_eQTA = Reduce(rbind, top_eQTA_list)

data_plot = data.frame(celltype = top_eQTA$celltype, num = top_eQTA$top_nearest, group = 'Top eGene is the nearest gene')
data_plot = rbind(data_plot, data.frame(celltype = top_eQTA$celltype, num = top_eQTA$top_not_nearest, group = 'Top eGene is not the nearest gene'))
data_plot$celltype = factor(data_plot$celltype, levels = levels)

p1 = ggplot(data_plot, aes(x = num, y = celltype, fill = group))+
  geom_col(position = 'fill', col = 'white') + 
  geom_text(aes(label = celltype, x = 0.01), hjust = 0, colour = "black")+
  theme_bw()+
  scale_fill_manual(values = c('Top eGene is the nearest gene' = 'dodgerblue3', 'Top eGene is not the nearest gene' = 'firebrick'))+
  scale_x_continuous(n.breaks = 10)+
  labs(y = 'Cell types', x = 'Proportions', fill = '')+
  theme(legend.position = 'bottom', axis.text.y = element_blank(), axis.ticks.y = element_blank())
pdf('../plots/eQTA/top_eQTA.pdf', width = 5, height = 5)
p1
dev.off()
############################
load('~/RWorkSpace/Duerr_bulk/eQTA/eQTA.RData')

bulk = me$cis$eqtls
bulk_sig = bulk[bulk$FDR < 0.05,]

check_overlap_bulk = function(celltype, ft1){
  ft2 = df_list[[celltype]]
  ft_shared = intersect(ft1, ft2)
  return(length(ft_shared)/length(ft2))
}

df_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  df <- check_peak(celltype)
  df_list[[celltype]] <- df
}

bulk_peak = sapply(levels[1:20], function(i){
  check_overlap_bulk(i,unique(bulk$snps))
})

df_list <- list()
for(i in 1:20){
  celltype <- levels[i]
  df <- check_gene(celltype)
  df_list[[celltype]] <- df
}

bulk_gene = sapply(levels[1:20], function(i){
  check_overlap_bulk(i,unique(bulk$gene))
})

pairs = p2g_celltype_link
for (i in 1:20){
  pairs[[i]]$pair = paste(pairs[[i]]$peak_chr_pos, pairs[[i]]$gene, sep = ":")
}

bulk$pair = paste(bulk$snps, bulk$gene, sep = ":")

pair_replicate = function(celltype){
  df = pairs[[celltype]]
  df_pair = merge(bulk, df, by = 'pair')
  table = data.frame(celltype = celltype, overlap = nrow(df_pair), dogma = nrow(df), bulk = nrow(bulk), replicate = as.numeric(table(df_pair$beta * df_pair$Correlation > 0)['TRUE']))
  return(table)
}

pair_rp_list = list()
for (i in 1:20){
  pair_rp_list[[levels[i]]] = pair_replicate(levels[i])
}
pair_rp = Reduce(rbind, pair_rp_list)
fwrite(pair_rp, file = '../plots/eQTA/bulk_replicate.txt', row.names = F, col.names = T)

################
load('../plots/eQTA/pseudobulk_eQTA_re.RData')

pseudobulk = me$cis$eqtls
pseudobulk$pair = paste(pseudobulk$snps, pseudobulk$gene, sep = ":")

length(intersect(unique(bulk$snps), unique(pseudobulk$snps)))#142537
length(intersect(unique(bulk$gene), unique(pseudobulk$gene)))#15594

df_pair = merge(bulk, pseudobulk, by = 'pair')
table = data.frame(overlap = nrow(df_pair), bulk = nrow(bulk), pseudobulk = nrow(pseudobulk), replicate = as.numeric(table(df_pair$beta.x * df_pair$beta.y > 0)['TRUE']))

################
read_pseudobulk = function(celltype){
  me = readRDS(paste0('../output/eQTA/pseudobulk_eQTA_', celltype, '.RDS'))
  tmp = me$cis$eqtls
  tmp$pair = paste(tmp$snps, tmp$gene, sep = ":")
  return(tmp)
}

pseudobulk = lapply(levels[1:20], function(i){read_pseudobulk(i)})
names(pseudobulk) = levels

pseudobulk_pair_replicate = function(celltype){
  df = pairs[[celltype]]
  bulk_df = pseudobulk[[celltype]]
  df_pair = merge(bulk_df, df, by = 'pair')
  table = data.frame(celltype = celltype, overlap = nrow(df_pair), dogma = nrow(df), pseudobulk = nrow(bulk_df), replicate = as.numeric(table(df_pair$beta * df_pair$Correlation > 0)['TRUE']),
  fdr_0.05 = as.numeric(table(bulk_df$FDR < 0.05)['TRUE']), fdr_0.01 = as.numeric(table(bulk_df$FDR < 0.01)['TRUE']), fdr_0.001 = as.numeric(table(bulk_df$FDR < 0.001)['TRUE']),
  fdr_0.05_rp = as.numeric(table(df_pair$FDR.x < 0.05 & df_pair$beta * df_pair$Correlation > 0)['TRUE']), fdr_0.01_rp = as.numeric(table(df_pair$FDR.x < 0.01 & df_pair$beta * df_pair$Correlation > 0)['TRUE']), fdr_0.001_rp = as.numeric(table(df_pair$FDR.x < 0.001 & df_pair$beta * df_pair$Correlation > 0)['TRUE']))
  return(table)
}

pseudobulk_pair_rp_list = lapply(levels[1:20], function(i){pseudobulk_pair_replicate(i)})
pseudobulk_pair_rp = Reduce(rbind, pseudobulk_pair_rp_list)

fwrite(pseudobulk_pair_rp, file = '../plots/eQTA/pseudobulk_replicate.txt', row.names = F, col.names = T)
