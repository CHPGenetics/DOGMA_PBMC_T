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

# snps <- read.csv('../output/SNP/SNP_disease.csv', row.names = 'X')
# 
# disease_snps <- snps
# 
# disease_snps$Disease_merged <- NA
# disease_snps$Group_merged <- NA
# for (i in 1:nrow(disease_snps)){
#   tmp <- disease_snps[disease_snps$chr_pos == disease_snps$chr_pos[i],]
#   disease_snps$Disease_merged[i] <- paste(tmp$Disease, collapse = ',')
#   disease_snps$Group_merged[i] <- paste(tmp$Group, collapse = ',')
# }
# 
# disease_snps <- disease_snps[!duplicated(disease_snps$chr_pos),]
# disease_snps <- disease_snps[,c('SNP', 'chr_38', 'pos_38', 'Disease_merged', 'Group_merged')]
# disease_snps_gr <- makeGRangesFromDataFrame(disease_snps, seqnames.field = 'chr_38', start.field = 'pos_38', end.field = 'pos_38', keep.extra.columns = T)
# 
# saveRDS(disease_snps_gr, file = '../output/SNP/disease_gr.RDS')
disease_snps_gr <- readRDS('../output/SNP/disease_gr.RDS')

# snps <- read.csv('../output/SNP/SNP_disease.csv', row.names = 'X')
# 
# disease_snps <- snps[snps$Group == 'Immune',]
# 
# disease_snps$Disease_merged <- NA
# disease_snps$Group_merged <- NA
# for (i in 1:nrow(disease_snps)){
#   tmp <- disease_snps[disease_snps$chr_pos == disease_snps$chr_pos[i],]
#   disease_snps$Disease_merged[i] <- paste(tmp$Disease, collapse = ',')
#   disease_snps$Group_merged[i] <- paste(tmp$Group, collapse = ',')
# }
# 
# disease_snps <- disease_snps[!duplicated(disease_snps$chr_pos),]
# disease_snps <- disease_snps[,c('SNP', 'chr_38', 'pos_38', 'Disease_merged', 'Group_merged')]
# disease_snps_gr <- makeGRangesFromDataFrame(disease_snps, seqnames.field = 'chr_38', start.field = 'pos_38', end.field = 'pos_38', keep.extra.columns = T)
# 
# saveRDS(disease_snps_gr, file = '../output/SNP/disease_immune_gr.RDS')
disease_snps_gr <- readRDS('../output/SNP/disease_immune_gr.RDS')

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

p2g_celltype_link <- list()
for(i in 1:9){
  celltype <- levels[i]
  name <- paste0('0',i,'_',celltype)
  p2g <- readRDS(paste0('../output/ArchR/DOGMA_filtered_multiome/', name, '/p2g_fdr_0.05_loop_replicate_peak_gene.RDS'))
  p2g_celltype_link[[celltype]] <- p2g$Peak2GeneLinks
}
for(i in 10:20){
  celltype <- levels[i]
  name <- paste0(i,'_',celltype)
  p2g <- readRDS(paste0('../output/ArchR/DOGMA_filtered_multiome/', name, '/p2g_fdr_0.05_loop_replicate_peak_gene.RDS'))
  p2g_celltype_link[[celltype]] <- p2g$Peak2GeneLinks
}

disease_celltype <- function(celltype, disease){
  p2g_tmp <- Reduce(c, p2g_celltype_link[celltype])
  overlap <- subsetByOverlaps(p2g_tmp, disease_snps_gr[str_detect(disease_snps_gr$Disease_merged, disease),])
  return(overlap)
}

meta = fread('../plots/SNP/CHEERS/enrich_disease_individual_meta_bonferroni.txt')

count_gene <- function(celltype, disease){
  gr <- disease_celltype(celltype, disease)
  head(table(gr$gene)[order(table(gr$gene),decreasing = T)], n = 20)
}

measure_gene <- function(celltype, disease, gene){
  gr <- disease_celltype(celltype, disease)
  gr1 <- gr[gr$gene == gene,]
  hist(gr1$distance, breaks = seq(-1e6, 1e6, 1e5))
}

count_gene_multi <- function(meta, disease){
  celltype <- str_split_fixed(meta[str_detect(meta$SNP, disease) & meta$p_val_adj < 0.05,]$name, '_', 2)[,2]
  celltype <- unique(celltype)
  print(celltype)
  gr <- disease_celltype(celltype, disease)
  head(table(gr$gene)[order(table(gr$gene),decreasing = T)], n = 20)
}

measure_gene_multi <- function(meta, disease, gene){
  celltype <- str_split_fixed(meta[str_detect(meta$SNP, disease) & meta$p_val_adj < 0.05,]$name, '_', 2)[,2]
  gr <- disease_celltype(celltype, disease)
  gr1 <- gr[gr$gene == gene,]
  hist(gr1$distance, breaks = seq(-1e6, 1e6, 1e5))
}

count_gene_multi(meta, 'Psoriasis')
#IFNLR1,RUNX3
measure_gene_multi(meta, 'Psoriasis', 'IFNLR1')
#0,8e5
measure_gene_multi(meta, 'Psoriasis', 'RUNX3')
#-8e5,1e5

count_gene_multi(meta, 'Behcets_disease')
#KLRC4
measure_gene_multi(meta, 'Behcets_disease', 'KLRC4')
#-1e5,0

count_gene_multi(meta, 'Alopecia_areata')
#IL21
measure_gene_multi(meta, 'Alopecia_areata', 'IL21')
#-1e5,1e5

count_gene_multi(meta, 'Ankylosing_spondylitis')
#GPR65
measure_gene_multi(meta, 'Ankylosing_spondylitis', 'GPR65')
#-1e5,1e5

count_gene_multi(meta, 'Systemic_lupus_erythematosus')
#BLK,CD44
measure_gene_multi(meta, 'Systemic_lupus_erythematosus', 'BLK')
#0,1e5
measure_gene_multi(meta, 'Systemic_lupus_erythematosus', 'CD44')
#-1e5,0

####################
# snp_all <- read.table('../output/SNP/SNP_IBD_all.txt', header = T)
# snp_lead <- read.table('../output/SNP/SNP_IBD_lead.txt', header = T)
# snp_p <- read.table('../output/SNP/SNP_IBD_p.txt', header = T)
# 
# disease_snps <- snp_all
# disease_snps$Type <- 'All'
# disease_snps$Type[disease_snps$SNP %in% snp_lead$SNP] <- paste(disease_snps$Type[disease_snps$SNP %in% snp_lead$SNP], 'Lead', sep = ',')
# disease_snps$Type[disease_snps$SNP %in% snp_p$SNP] <- paste(disease_snps$Type[disease_snps$SNP %in% snp_p$SNP], 'Causal', sep = ',')
# 
# disease_snps <- disease_snps[,c('SNP', 'chr_38', 'pos_38', 'Type', 'trait_uc', 'trait_cd')]
# disease_snps$Disease_merged <- disease_snps$Type
# disease_snps_gr <- makeGRangesFromDataFrame(disease_snps, seqnames.field = 'chr_38', start.field = 'pos_38', end.field = 'pos_38', keep.extra.columns = T)
# 
# saveRDS(disease_snps_gr, file = '../output/SNP/ibd_gr.RDS')
disease_snps_gr <- readRDS('../output/SNP/ibd_gr.RDS')
#disease_snps_gr <- disease_snps_gr[disease_snps_gr$Type %in% c('All,Lead','All,Lead,Causal')]
meta_ibd = fread('../plots/SNP/CHEERS/enrich_disease_IBD_individual_meta_bonferroni.txt')

count_gene(levels[6] , 'All')
measure_gene(levels[6] , 'All', 'PTGER4')
#-2e5,0
measure_gene(levels[6] , 'All', 'CREM')
#-2e5,1e5

count_gene_multi(meta_ibd, 'All')
#PTGER4, LNPEP, TNFRSF14, IL2RA
count_gene_multi(meta_ibd, 'Lead')
#CCL20, CXCR5, ITGAL

measure_gene_multi(meta_ibd, 'All', 'PTGER4')
#-5e5,0
measure_gene_multi(meta_ibd, 'All', 'LNPEP')
#-1e5,1e5
measure_gene_multi(meta_ibd, 'All', 'TNFRSF14')
#-1e5,1e5
measure_gene_multi(meta_ibd, 'All', 'IL2RA')
#-1e5,0

measure_gene_multi(meta_ibd, 'Lead', 'CCL20')
#-1e5,0
measure_gene_multi(meta_ibd, 'Lead', 'CXCR5')
#-1e5,1e5
measure_gene_multi(meta_ibd, 'Lead', 'ITGAL')
#-1e5,0

meta_ibd_ibd <- meta_ibd[meta_ibd$group == 'IBD',]
count_gene_multi(meta_ibd_ibd, 'All')
#PTGER4, LNPEP, TNFRSF14, IL2RA
count_gene_multi(meta_ibd_ibd, 'Lead')
#CCL20, CXCR5, ITGAL

meta_ibd_cd <- meta_ibd[meta_ibd$group == 'CD',]
count_gene_multi(meta_ibd_cd, 'All')
#PTGER4, LNPEP, TNFRSF14, IL2RA
count_gene_multi(meta_ibd_cd, 'Lead')
#CCL20, CXCR5, ITGAL

meta_ibd_uc <- meta_ibd[meta_ibd$group == 'UC',]
count_gene_multi(meta_ibd_uc, 'All')
#PTGER4, LNPEP, TNFRSF14, IL2RA
count_gene_multi(meta_ibd_uc, 'Lead')
#CCL20, CXCR5, ITGAL