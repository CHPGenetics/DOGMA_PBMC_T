library(data.table)
library(stringr)
library(liftOver)
library(SNPlocs.Hsapiens.dbSNP141.GRCh38)
library(rtracklayer)
library(Signac)
library(cowplot)
library(ggsci)
library(ggpubr)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

p2g <- readRDS('../output/ArchR/DOGMA_filtered_multiome/p2g_fdr_0.05.RDS')

snps <- read.csv('../output/SNP/SNP_disease.csv', row.names = 'X')

disease_snps <- snps

disease_snps$Disease_merged <- NA
for (i in 1:nrow(disease_snps)){
  tmp <- disease_snps[disease_snps$chr_pos == disease_snps$chr_pos[i],]
  disease_snps$Disease_merged[i] <- paste(tmp$Disease, collapse = ',')
}

disease_snps <- disease_snps[!duplicated(disease_snps$chr_pos),]
disease_snps <- disease_snps[,c('SNP', 'chr_38', 'pos_38', 'Disease_merged')]
disease_snps <- makeGRangesFromDataFrame(disease_snps, seqnames.field = 'chr_38', start.field = 'pos_38', end.field = 'pos_38', keep.extra.columns = T)

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

peak_df_sub <- peak_df[!duplicated(peak_df$chr_pos),]
rownames(peak_df_sub) <- peak_df_sub$chr_pos
peak_df_sub_gr <- makeGRangesFromDataFrame(peak_df_sub, seqnames.field = 'chr', start.field = 'start', end.field = 'end')

overlap_idx <- findOverlaps(peak_df_sub_gr, disease_snps)
overlap_peak <- peak_df_sub_gr[overlap_idx@from]
overlap_snp <- disease_snps[overlap_idx@to]

overlap_peak_df <- data.frame(chr = seqnames(overlap_peak), start = start(overlap_peak), end = end(overlap_peak))
overlap_peak_df$chr_pos <- paste(overlap_peak_df$chr, overlap_peak_df$start, overlap_peak_df$end, sep = '-')

overlap_snp_df <- data.frame(chr = seqnames(overlap_snp), pos = start(overlap_snp), SNP = overlap_snp$SNP, Disease_merged = overlap_snp$Disease_merged)
overlap_snp_df$chr_pos <- paste(overlap_snp_df$chr, overlap_snp_df$pos, sep = '-')

overlap_df <- cbind(overlap_peak_df, overlap_snp_df)
colnames(overlap_df)[1:4] <- paste0('peak_', colnames(overlap_df)[1:4])
colnames(overlap_df)[c(5,6,ncol(overlap_df))] <- paste0('snp_', colnames(overlap_df)[c(5,6,ncol(overlap_df))])
overlap_df <- overlap_df[,-c(1:3)]

overlap_df$snp_chr_merged <- NA
overlap_df$snp_pos_merged <- NA
overlap_df$SNP_merged <- NA
overlap_df$Disease_merged_merged <- NA
overlap_df$snp_chr_pos_merged <- NA
for (i in 1:nrow(overlap_df)){
  tmp <- overlap_df[overlap_df$peak_chr_pos == overlap_df$peak_chr_pos[i],]
  overlap_df$snp_chr_merged[i] <- paste(tmp$snp_chr, collapse = '|')
  overlap_df$snp_pos_merged[i] <- paste(tmp$snp_pos, collapse = '|')
  overlap_df$SNP_merged[i] <- paste(tmp$SNP, collapse = '|')
  overlap_df$Disease_merged_merged[i] <- paste(tmp$Disease_merged, collapse = '|')
  overlap_df$snp_chr_pos_merged[i] <- paste(tmp$snp_chr, collapse = '|')
}

p2g_disease_peak <- as.data.frame(str_split_fixed(overlap_df$peak_chr_pos, '-', 3))
write.table(p2g_disease_peak, file = '../output/SNP/p2g_disease_peak.bed', quote = F, row.names = F, col.names = F)
p2g_disease_snp <- overlap_df[,c('snp_chr', 'snp_pos', 'snp_pos', 'SNP')]
write.table(p2g_disease_snp, file = '../output/SNP/p2g_disease_snp.bed', quote = F, row.names = F, col.names = F)

overlap_df <- overlap_df[!duplicated(overlap_df$peak_chr_pos),]
overlap_df <- overlap_df[,-c(2:6)]
p2g_disease <- merge(p2g_df, overlap_df, by = 'peak_chr_pos', all.x = T)

write.table(p2g_disease, file = '../output/SNP/p2g_disease.txt', quote = F, row.names = F, col.names = T)

p2g_disease$peak_pos <- (p2g_disease$peak_start + p2g_disease$peak_end)/2
p2g_disease$min <- pmin(p2g_disease$peak_pos, p2g_disease$gene_pos)
p2g_disease$max <- pmax(p2g_disease$peak_pos, p2g_disease$gene_pos)
p2g_disease$score <- 0.5
p2g_disease$score[!is.na(p2g_disease$SNP_merged)] <- 1
p2g_disease$group <- 2*p2g_disease$score

links <- makeGRangesFromDataFrame(p2g_disease, seqnames.field = 'peak_chr', start.field = 'min', end.field = 'max', keep.extra.columns = T)
saveRDS(links, file = '../output/SNP/links_disease.RDS')
#################
snp_all <- read.table('../output/SNP/SNP_IBD_all.txt', header = T)
snp_lead <- read.table('../output/SNP/SNP_IBD_lead.txt', header = T)
snp_p <- read.table('../output/SNP/SNP_IBD_p.txt', header = T)

disease_snps <- snp_all
disease_snps$Type <- 'All'
disease_snps$Type[disease_snps$SNP %in% snp_lead$SNP] <- paste(disease_snps$Type[disease_snps$SNP %in% snp_lead$SNP], 'Lead', sep = ',')
disease_snps$Type[disease_snps$SNP %in% snp_p$SNP] <- paste(disease_snps$Type[disease_snps$SNP %in% snp_p$SNP], 'Causal', sep = ',')

disease_snps <- disease_snps[,c('SNP', 'chr_38', 'pos_38', 'Type', 'trait_uc', 'trait_cd')]
disease_snps <- makeGRangesFromDataFrame(disease_snps, seqnames.field = 'chr_38', start.field = 'pos_38', end.field = 'pos_38', keep.extra.columns = T)

# peak <- metadata(p2g)$peakSet
# gene <- metadata(p2g)$geneSet
# p2g_df <- as.data.frame(p2g)
# 
# peak_sub <- peak[p2g_df$idxATAC]
# gene_sub <- gene[p2g_df$idxRNA]
# peak_df <- data.frame(chr = seqnames(peak_sub), start = start(peak_sub), end = end(peak_sub))
# gene_df <- data.frame(chr = seqnames(gene_sub), pos = start(gene_sub), gene = gene_sub$name)
# peak_df$chr_pos <- paste(peak_df$chr, peak_df$start, peak_df$end, sep = '-')
# 
# p2g_df$peak_chr <- peak_df$chr
# p2g_df$peak_start <- peak_df$start
# p2g_df$peak_end <- peak_df$end
# p2g_df$peak_chr_pos <- peak_df$chr_pos
# p2g_df$gene_chr <- gene_df$chr
# p2g_df$gene_pos <- gene_df$pos
# p2g_df$gene <- gene_df$gene
# length(unique(p2g_df$peak_chr_pos))#26874
# length(unique(p2g_df$gene))#6652
# 
# peak_df_sub <- peak_df[!duplicated(peak_df$chr_pos),]
# rownames(peak_df_sub) <- peak_df_sub$chr_pos
# peak_df_sub_gr <- makeGRangesFromDataFrame(peak_df_sub, seqnames.field = 'chr', start.field = 'start', end.field = 'end')

overlap_idx <- findOverlaps(peak_df_sub_gr, disease_snps)
overlap_peak <- peak_df_sub_gr[overlap_idx@from]
overlap_snp <- disease_snps[overlap_idx@to]

overlap_peak_df <- data.frame(chr = seqnames(overlap_peak), start = start(overlap_peak), end = end(overlap_peak))
overlap_peak_df$chr_pos <- paste(overlap_peak_df$chr, overlap_peak_df$start, overlap_peak_df$end, sep = '-')

overlap_snp_df <- data.frame(chr = seqnames(overlap_snp), pos = start(overlap_snp), SNP = overlap_snp$SNP, Type = overlap_snp$Type, trait_uc = overlap_snp$trait_uc, trait_cd = overlap_snp$trait_cd)
overlap_snp_df$chr_pos <- paste(overlap_snp_df$chr, overlap_snp_df$pos, sep = '-')

overlap_df <- cbind(overlap_peak_df, overlap_snp_df)
colnames(overlap_df)[1:4] <- paste0('peak_', colnames(overlap_df)[1:4])
colnames(overlap_df)[c(5,6,ncol(overlap_df))] <- paste0('snp_', colnames(overlap_df)[c(5,6,ncol(overlap_df))])
overlap_df <- overlap_df[,-c(1:3)]

overlap_df$snp_chr_merged <- NA
overlap_df$snp_pos_merged <- NA
overlap_df$SNP_merged <- NA
overlap_df$Type_merged <- NA
overlap_df$trait_uc_merged <- NA
overlap_df$trait_cd_merged <- NA
overlap_df$snp_chr_pos_merged <- NA
for (i in 1:nrow(overlap_df)){
  tmp <- overlap_df[overlap_df$peak_chr_pos == overlap_df$peak_chr_pos[i],]
  overlap_df$snp_chr_merged[i] <- paste(tmp$snp_chr, collapse = '|')
  overlap_df$snp_pos_merged[i] <- paste(tmp$snp_pos, collapse = '|')
  overlap_df$SNP_merged[i] <- paste(tmp$SNP, collapse = '|')
  overlap_df$Type_merged[i] <- paste(tmp$Type, collapse = '|')
  overlap_df$trait_uc_merged[i] <- paste(tmp$trait_uc, collapse = '|')
  overlap_df$trait_cd_merged[i] <- paste(tmp$trait_cd, collapse = '|')
  overlap_df$snp_chr_pos_merged[i] <- paste(tmp$snp_chr, collapse = '|')
}

p2g_ibd_peak <- as.data.frame(str_split_fixed(overlap_df$peak_chr_pos, '-', 3))
write.table(p2g_ibd_peak, file = '../output/SNP/p2g_ibd_peak.bed', quote = F, row.names = F, col.names = F)
p2g_ibd_snp <- overlap_df[,c('snp_chr', 'snp_pos', 'snp_pos', 'SNP')]
write.table(p2g_ibd_snp, file = '../output/SNP/p2g_ibd_snp.bed', quote = F, row.names = F, col.names = F)

overlap_df <- overlap_df[!duplicated(overlap_df$peak_chr_pos),]
overlap_df <- overlap_df[,-c(2:8)]
p2g_ibd <- merge(p2g_df, overlap_df, by = 'peak_chr_pos', all.x = T)

write.table(p2g_ibd, file = '../output/SNP/p2g_ibd.txt', quote = F, row.names = F, col.names = T)

p2g_ibd$peak_pos <- (p2g_ibd$peak_start + p2g_ibd$peak_end)/2
p2g_ibd$min <- pmin(p2g_ibd$peak_pos, p2g_ibd$gene_pos)
p2g_ibd$max <- pmax(p2g_ibd$peak_pos, p2g_ibd$gene_pos)
p2g_ibd$score <- 0.5
p2g_ibd$score[!is.na(p2g_ibd$SNP_merged)] <- 1
p2g_ibd$group <- 2*p2g_ibd$score

links <- makeGRangesFromDataFrame(p2g_ibd, seqnames.field = 'peak_chr', start.field = 'min', end.field = 'max', keep.extra.columns = T)
saveRDS(links, file = '../output/SNP/links_ibd.RDS')
#################
snp_all <- read.table('../output/SNP/SNP_IBD_all.txt', header = T)
snp_lead <- read.table('../output/SNP/SNP_IBD_lead.txt', header = T)
snp_p <- read.table('../output/SNP/SNP_IBD_p.txt', header = T)

disease_snps <- snp_all
disease_snps$Type <- 'All'
disease_snps$Type[disease_snps$SNP %in% snp_lead$SNP] <- paste(disease_snps$Type[disease_snps$SNP %in% snp_lead$SNP], 'Lead', sep = ',')
disease_snps$Type[disease_snps$SNP %in% snp_p$SNP] <- paste(disease_snps$Type[disease_snps$SNP %in% snp_p$SNP], 'Causal', sep = ',')

disease_snps <- disease_snps[,c('SNP', 'chr_38', 'pos_38', 'Type', 'trait_uc', 'trait_cd')]
disease_snps <- makeGRangesFromDataFrame(disease_snps, seqnames.field = 'chr_38', start.field = 'pos_38', end.field = 'pos_38', keep.extra.columns = T)

load('../output/ArchR/DOGMA_filtered_multiome/PeakMatrix.RData')

peaks_mat_df <- data.frame(chr=seqnames(peaks_mat),
                           start=start(peaks_mat),
                           end=end(peaks_mat),
                           distToTSS = rowRanges(peaks_mat)$distToTSS,
                           GC = rowRanges(peaks_mat)$GC)
rownames(peaks_mat_df) <- GRangesToString(rowRanges(peaks_mat), sep = c("-", "-"))

peak_df_sub_gr <- makeGRangesFromDataFrame(peaks_mat_df, seqnames.field = 'chr', start.field = 'start', end.field = 'end')

overlap_idx <- findOverlaps(peak_df_sub_gr, disease_snps)
overlap_peak <- peak_df_sub_gr[overlap_idx@from]
overlap_snp <- disease_snps[overlap_idx@to]

overlap_peak_df <- data.frame(chr = seqnames(overlap_peak), start = start(overlap_peak), end = end(overlap_peak))
overlap_peak_df$chr_pos <- paste(overlap_peak_df$chr, overlap_peak_df$start, overlap_peak_df$end, sep = '-')

overlap_snp_df <- data.frame(chr = seqnames(overlap_snp), pos = start(overlap_snp), SNP = overlap_snp$SNP, Type = overlap_snp$Type, trait_uc = overlap_snp$trait_uc, trait_cd = overlap_snp$trait_cd)
overlap_snp_df$chr_pos <- paste(overlap_snp_df$chr, overlap_snp_df$pos, sep = '-')

overlap_df <- cbind(overlap_peak_df, overlap_snp_df)
colnames(overlap_df)[1:4] <- paste0('peak_', colnames(overlap_df)[1:4])
colnames(overlap_df)[c(5,6,ncol(overlap_df))] <- paste0('snp_', colnames(overlap_df)[c(5,6,ncol(overlap_df))])
overlap_df <- overlap_df[,-c(1:3)]

overlap_df$snp_chr_merged <- NA
overlap_df$snp_pos_merged <- NA
overlap_df$SNP_merged <- NA
overlap_df$Type_merged <- NA
overlap_df$trait_uc_merged <- NA
overlap_df$trait_cd_merged <- NA
overlap_df$snp_chr_pos_merged <- NA
for (i in 1:nrow(overlap_df)){
  tmp <- overlap_df[overlap_df$peak_chr_pos == overlap_df$peak_chr_pos[i],]
  overlap_df$snp_chr_merged[i] <- paste(tmp$snp_chr, collapse = '|')
  overlap_df$snp_pos_merged[i] <- paste(tmp$snp_pos, collapse = '|')
  overlap_df$SNP_merged[i] <- paste(tmp$SNP, collapse = '|')
  overlap_df$Type_merged[i] <- paste(tmp$Type, collapse = '|')
  overlap_df$trait_uc_merged[i] <- paste(tmp$trait_uc, collapse = '|')
  overlap_df$trait_cd_merged[i] <- paste(tmp$trait_cd, collapse = '|')
  overlap_df$snp_chr_pos_merged[i] <- paste(tmp$snp_chr_pos, collapse = '|')
}

p2g_ibd_peak <- as.data.frame(str_split_fixed(overlap_df$peak_chr_pos, '-', 3))
write.table(p2g_ibd_peak, file = '../output/SNP/overlap_ibd_peak.bed', quote = F, row.names = F, col.names = F)
p2g_ibd_snp <- overlap_df[,c('snp_chr', 'snp_pos', 'snp_pos', 'SNP')]
write.table(p2g_ibd_snp, file = '../output/SNP/overlap_ibd_snp.bed', quote = F, row.names = F, col.names = F)

###############
peak_annotation <- readRDS('../output/pseudo_bulk/peak_annotation.RDS')

overlap_df_annotation <- peak_annotation[overlap_df$peak_chr_pos,]
overlap_df <- cbind(overlap_df, overlap_df_annotation)
overlap_df_dedup <- overlap_df[!duplicated(overlap_df$peak_chr_pos),]

write.table(overlap_df_dedup, file = '../output/SNP/overlap_ibd_peak_snp.txt', quote = F, row.names = T, col.names = T, sep = '\t')

table(overlap_df$gene_name)[order(table(overlap_df$gene_name), decreasing = T)]
table(overlap_df_dedup$gene_name)[order(table(overlap_df_dedup$gene_name), decreasing = T)]
