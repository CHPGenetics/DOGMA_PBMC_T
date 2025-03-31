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
  
  #p2g_df <- p2g_df[,-c(1:2)]
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

pairs = p2g_celltype_link
for (i in 1:20){
  pairs[[i]]$pair = paste(pairs[[i]]$peak_chr_pos, pairs[[i]]$gene, sep = ":")
}
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
  df_pair_rp = df_pair[df_pair$beta * df_pair$Correlation > 0,]
  df_pair_rp_1 = df_pair_rp[,-c(1:7)]
  colnames(df_pair_rp_1)[c(4,13)] = c('FDR', 'gene')
  df_pair_rp_1$pair = paste(df_pair_rp_1$peak_chr_pos, df_pair_rp_1$gene, sep = ":")
  return(df_pair_rp_1)
}

pairs_replicate = lapply(levels[1:20], function(i){pseudobulk_pair_replicate(i)})
names(pairs_replicate) = levels

###############
check_pair = function(celltype){
  df <- pairs_replicate[[celltype]]
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
  geom_text(aes(label = celltype, x = 1e3), hjust = 0, colour = "black")+
  theme_bw()+
  scale_fill_manual(values = c('Positive' = 'firebrick', 'Negative' = 'dodgerblue3'))+
  scale_x_continuous(n.breaks = 12)+
  labs(y = 'Cell types', x = 'Peak-to-gene linkages within 1M bps', fill = '')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'bottom', axis.text.y = element_blank(), axis.ticks.y = element_blank())
pdf('../plots/eQTA_replicate/pairs.pdf', width = 3, height = 5)
p1
dev.off()

df$celltype = factor(df$celltype, levels = levels)
p1 = ggplot(df, aes(x = genes, y = celltype, fill = celltype))+
  geom_col(position = 'stack', col = 'white')+
  geom_text(aes(label = celltype, x = 200), hjust = 0, colour = "black")+
  theme_bw()+
  scale_x_continuous(n.breaks = 10)+
  labs(y = 'Cell types', x = '# of eGenes', fill = '')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', axis.text.y = element_blank(), axis.ticks.y = element_blank())
pdf('../plots/eQTA_replicate/genes.pdf', width = 3, height = 5)
p1
dev.off()

p1 = ggplot(df, aes(x = peaks, y = celltype, fill = celltype))+
  geom_col(position = 'stack', col = 'white')+
  geom_text(aes(label = celltype, x = 1000), hjust = 0, colour = "black")+
  theme_bw()+
  scale_x_continuous(n.breaks = 7)+
  labs(y = 'Cell types', x = '# of ePeaks', fill = '')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', axis.text.y = element_blank(), axis.ticks.y = element_blank())
pdf('../plots/eQTA_replicate/peaks.pdf', width = 3, height = 5)
p1
dev.off()

###############
check_gene = function(celltype){
  df <- pairs_replicate[[celltype]]
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
pdf('../plots/eQTA_replicate/shared_gene.pdf', width = 8, height = 5)
Heatmap(mat*100, cluster_rows = F, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = "%"), col = viridis(100, option = "C", direction = -1))
dev.off()

###############
check_peak = function(celltype){
  df <- pairs_replicate[[celltype]]
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

pdf('../plots/eQTA_replicate/shared_peak.pdf', width = 8, height = 5)
Heatmap(mat*100, cluster_rows = F, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = "%"), col = viridis(100, option = "C", direction = -1))
dev.off()

###############
#peak_annotation <- readRDS('../output/pseudo_bulk/peak_annotation.RDS')
#peak_annotation$name = rownames(peak_annotation)

anno <- fread('../output/pseudo_bulk/peaks.annotation_sorted.txt')[,-1]
anno$name <- paste(anno$Chr, anno$Start-1, anno$End, sep = '-')

check_top_eQTA = function(celltype){
  tmp = pairs_replicate[[celltype]]
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
data_plot$group = factor(data_plot$group, levels = c('Top eGene is the nearest gene', 'Top eGene is not the nearest gene'))

p1 = ggplot(data_plot, aes(x = num, y = celltype, fill = group))+
  geom_col(position = 'fill', col = 'white') + 
  geom_text(aes(label = celltype, x = 0.01), hjust = 0, colour = "black")+
  theme_bw()+
  scale_fill_manual(values = c('Top eGene is not the nearest gene' = 'firebrick', 'Top eGene is the nearest gene' = 'dodgerblue3'))+
  scale_x_continuous(n.breaks = 10)+
  labs(y = 'Cell types', x = 'Proportions', fill = '')+
  theme(legend.position = 'bottom', axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  guides(fill = guide_legend(ncol = 1, override.aes = list(size=4)))
pdf('../plots/eQTA_replicate/top_eQTA.pdf', width = 3, height = 5)
p1
dev.off()

####################
pseudobulk_pair_replicate = function(celltype){
  df = pairs[[celltype]]
  bulk_df = pseudobulk[[celltype]]
  df_pair = merge(bulk_df, df, by = 'pair')
  df_pair_rp = df_pair[df_pair$beta * df_pair$Correlation > 0,]
  df_pair_rp_1 = df_pair_rp[,-c(1:7)]
  colnames(df_pair_rp_1)[c(4,13)] = c('FDR', 'gene')
  df_pair_rp_1$pair = paste(df_pair_rp_1$peak_chr_pos, df_pair_rp_1$gene, sep = ":")
  return(df_pair_rp_1)
}

pairs_replicate = lapply(levels[1:20], function(i){pseudobulk_pair_replicate(i)})
names(pairs_replicate) = levels

form_loop = function(celltype){
    df = pairs_replicate[[celltype]]
    df_loop = data.frame(chr = df$gene_chr, gene = df$gene, gene_pos = df$gene_pos, peak_pos = (df$peak_start + df$peak_end)/2, value = df$Correlation, FDR = df$FDR)
    df_loop$start = pmin(df_loop$gene_pos, df_loop$peak_pos)
    df_loop$end = pmax(df_loop$gene_pos, df_loop$peak_pos)
    df_loop_re = df_loop
    grange = makeGRangesFromDataFrame(df_loop_re, seqnames.field = 'chr', start.field = 'start', end.field = 'end', keep.extra.columns = T)
    loop = list()
    loop[['Peak2GeneLinks']] = grange
    return(loop)
}

for (i in 1:20){
    celltype = levels[i]
    print(celltype)
    p2g = form_loop(celltype)
    if (i < 10){
        name <- paste0('0',i,'_',celltype)
    }
    else if (i >= 10) {
       name <- paste0(i,'_',celltype)
    }
    saveRDS(p2g, paste0('../output/ArchR/DOGMA_filtered_multiome/', name, '/p2g_fdr_0.001_loop_replicate.RDS'))
}

form_loop_peak_gene = function(list, celltype){
  df = list[[celltype]]
  df_loop = data.frame(chr = df$gene_chr, value = df$Correlation, FDR = df$FDR, peak = df$peak_chr_pos, gene = df$gene, peak_pos = (df$peak_start + df$peak_end)/2, gene_pos = df$gene_pos, start = df$peak_start, end = df$peak_end)
  df_loop$distance = df_loop$peak_pos - df_loop$gene_pos
  df_loop_re = df_loop
  grange = makeGRangesFromDataFrame(df_loop_re, seqnames.field = 'chr', start.field = 'start', end.field = 'end', keep.extra.columns = T)
  loop = list()
  loop[['Peak2GeneLinks']] = grange
  return(loop)
}

for (i in 1:20){
  celltype = levels[i]
  print(celltype)
  p2g = form_loop_peak_gene(pairs_replicate, celltype)
  if (i < 10){
    name <- paste0('0',i,'_',celltype)
  }
  else if (i >= 10) {
    name <- paste0(i,'_',celltype)
  }
  saveRDS(p2g, paste0('../output/ArchR/DOGMA_filtered_multiome/', name, '/p2g_fdr_0.001_loop_replicate_peak_gene.RDS'))
}

#################
p2g_celltype <- lapply(levels[1:20], function(i){
  gr <- form_loop_peak_gene(p2g_celltype_link,i)[['Peak2GeneLinks']]
  return(gr)
})
names(p2g_celltype) <- levels

p2g_celltype_replicate <- lapply(levels[1:20], function(i){
  gr <- form_loop_peak_gene(pairs_replicate,i)[['Peak2GeneLinks']]
  return(gr)
})
names(p2g_celltype_replicate) <- levels

ibd_gr <- readRDS('../output/SNP/ibd_gr.RDS')
disease_gr <- readRDS('../output/SNP/disease_gr.RDS')

disease_celltype <- function(list, celltype, snp, type, disease){
  overlap <- subsetByOverlaps(list[[celltype]], snp[str_detect(snp@elementMetadata[,type], disease),])
  return(overlap)
}

count <- sapply(levels[1:20], function(i){
  length(disease_celltype(p2g_celltype, i, disease_gr, 'Group_merged', 'Immune'))
})

count_replicate <- sapply(levels[1:20], function(i){
  length(disease_celltype(p2g_celltype_replicate, i, disease_gr, 'Group_merged', 'Immune'))
})

View(cbind(count, count_replicate))

count_ibd <- sapply(levels[1:20], function(i){
  length(disease_celltype(p2g_celltype, i, ibd_gr, 'Disease_merged', ''))
})

count_replicate_ibd <- sapply(levels[1:20], function(i){
  length(disease_celltype(p2g_celltype_replicate, i, ibd_gr, 'Disease_merged', ''))
})

View(cbind(count_ibd, count_replicate_ibd))
