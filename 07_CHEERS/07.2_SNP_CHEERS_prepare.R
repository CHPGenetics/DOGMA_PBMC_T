library(data.table)
library(stringr)
library(liftOver)
library(SNPlocs.Hsapiens.dbSNP141.GRCh38)
library(rtracklayer)
library(Signac)
library(cowplot)
library(ggsci)
library(ggpubr)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')

snps <- read.csv('../output/SNP/SNP_disease.csv', row.names = 'X')

for(i in names(table(snps$Disease))){
  disease_snps <- snps[snps$Disease == i,]
  disease_snps <- disease_snps[!duplicated(disease_snps$chr_pos),]
  snps_out <- disease_snps[,c(45,1,46)]
  colnames(snps_out) <- c('Chrom', 'SNP', 'BP')
  dir.create(file.path('../output/SNP/CHEERS/Disease_SNP/', i))
  write.table(snps_out, file = paste0('../output/SNP/CHEERS/Disease_SNP/', i, '/SNP.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
  
  snps_out <- disease_snps[,c(1,45,46)]
  colnames(snps_out) <- c('name', 'chr', 'pos')
  write.table(snps_out, file = paste0('../output/SNP/CHEERS/Disease_SNP/', i, '/SNP_list.txt'), quote = F, row.names = F, col.names = F, sep = "\t")
}

for(i in names(table(snps$Group_upper))){
  disease_snps <- snps[snps$Group_upper == i,]
  disease_snps <- disease_snps[!duplicated(disease_snps$chr_pos),]
  snps_out <- disease_snps[,c(45,1,46)]
  colnames(snps_out) <- c('Chrom', 'SNP', 'BP')
  dir.create(file.path('../output/SNP/CHEERS/Disease_SNP/', i))
  write.table(snps_out, file = paste0('../output/SNP/CHEERS/Disease_SNP/', i, '/SNP.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
  
  snps_out <- disease_snps[,c(1,45,46)]
  colnames(snps_out) <- c('name', 'chr', 'pos')
  write.table(snps_out, file = paste0('../output/SNP/CHEERS/Disease_SNP/', i, '/SNP_list.txt'), quote = F, row.names = F, col.names = F, sep = "\t")
}

All <- read.table('../output/SNP/SNP_IBD_all.txt', header = T)
Lead <- read.table('../output/SNP/SNP_IBD_lead.txt', header = T)
Causal <- read.table('../output/SNP/SNP_IBD_p.txt', header = T)

for(j in c('All', 'Lead', 'Causal')){
  disease_snps <- get(j)
  snps_out <- disease_snps[,c(ncol(disease_snps)-2,1,ncol(disease_snps)-1)]
  colnames(snps_out) <- c('Chrom', 'SNP', 'BP')
  
  dir.create(file.path('../output/SNP/CHEERS/IBD_SNP/', paste0(j, '_IBD')))
  write.table(snps_out, file = paste0('../output/SNP/CHEERS/IBD_SNP/', paste0(j, '_IBD'), '/SNP.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
  
  snps_out <- disease_snps[,c(1,ncol(disease_snps)-2,ncol(disease_snps)-1)]
  colnames(snps_out) <- c('name', 'chr', 'pos')
  write.table(snps_out, file = paste0('../output/SNP/CHEERS/IBD_SNP/', paste0(j, '_IBD'), '/SNP_list.txt'), quote = F, row.names = F, col.names = F, sep = "\t")
  
  disease_snps <- get(j)
  disease_snps <- disease_snps[disease_snps$trait_uc == 1,]
  snps_out <- disease_snps[,c(ncol(disease_snps)-2,1,ncol(disease_snps)-1)]
  colnames(snps_out) <- c('Chrom', 'SNP', 'BP')
  
  dir.create(file.path('../output/SNP/CHEERS/IBD_SNP/', paste0(j, '_UC')))
  write.table(snps_out, file = paste0('../output/SNP/CHEERS/IBD_SNP/', paste0(j, '_UC'), '/SNP.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
  
  snps_out <- disease_snps[,c(1,ncol(disease_snps)-2,ncol(disease_snps)-1)]
  colnames(snps_out) <- c('name', 'chr', 'pos')
  write.table(snps_out, file = paste0('../output/SNP/CHEERS/IBD_SNP/', paste0(j, '_UC'), '/SNP_list.txt'), quote = F, row.names = F, col.names = F, sep = "\t")
  
  disease_snps <- get(j)
  disease_snps <- disease_snps[disease_snps$trait_cd == 1,]
  snps_out <- disease_snps[,c(ncol(disease_snps)-2,1,ncol(disease_snps)-1)]
  colnames(snps_out) <- c('Chrom', 'SNP', 'BP')
  
  dir.create(file.path('../output/SNP/CHEERS/IBD_SNP/', paste0(j, '_CD')))
  write.table(snps_out, file = paste0('../output/SNP/CHEERS/IBD_SNP/', paste0(j, '_CD'), '/SNP.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
  
  snps_out <- disease_snps[,c(1,ncol(disease_snps)-2,ncol(disease_snps)-1)]
  colnames(snps_out) <- c('name', 'chr', 'pos')
  write.table(snps_out, file = paste0('../output/SNP/CHEERS/IBD_SNP/', paste0(j, '_CD'), '/SNP_list.txt'), quote = F, row.names = F, col.names = F, sep = "\t")
}

load('../output/ArchR/DOGMA_filtered_multiome/PeakMatrix.RData')
markerList <- readRDS('../output/ArchR/DOGMA_filtered_multiome/markerList.RDS')
markerList <- markerList[-9]

peaks_mat_df <- data.frame(chr=seqnames(peaks_mat),
                           start=start(peaks_mat),
                           end=end(peaks_mat),
                           distToTSS = rowRanges(peaks_mat)$distToTSS,
                           GC = rowRanges(peaks_mat)$GC)
rownames(peaks_mat_df) <- GRangesToString(rowRanges(peaks_mat), sep = c("-", "-"))

meta <- as.data.frame(colData(peaks_mat))
mat <- assay(peaks_mat)

rownames(mat) <- rownames(peaks_mat_df)

select_peaks <- function(name){
  meta_sub <- meta[meta$celltype_id == name,]
  mat_sub <- mat[,rownames(meta_sub)]
  mat_sub_count <- rowSums(mat_sub)
  peak_out <- peaks_mat_df[,c(1:3)]
  peak_out$count <- mat_sub_count
  
  dir.create(file.path('../output/SNP/CHEERS/AllPeaks/', 'Merged'))
  dir.create(file.path('../output/SNP/CHEERS/AllPeaks/', 'Individual'))
  fwrite(peak_out, file = paste0('../output/SNP/CHEERS/AllPeaks/Merged/', name, '_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
  
  for(i in names(table(meta$sample))){
    meta_sub_sub <- meta_sub[meta_sub$sample == i,]
    mat_sub_sub <- mat_sub[,rownames(meta_sub_sub)]
    mat_sub_sub_count <- rowSums(mat_sub_sub)
    peak_out <- peaks_mat_df[,c(1:3)]
    peak_out$count <- mat_sub_sub_count
    
    fwrite(peak_out, file = paste0('../output/SNP/CHEERS/AllPeaks/Individual/', name, '_', i,'_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
  }
}

for(j in names(table(meta$celltype_id))){
  select_peaks(j)
}

select_peaks_condition <- function(name){
  meta_sub <- meta[meta$celltype_id == name,]
  mat_sub <- mat[,rownames(meta_sub)]
  # mat_sub_count <- rowSums(mat_sub)
  # peak_out <- peaks_mat_df[,c(1:3)]
  # peak_out$count <- mat_sub_count
  
  dir.create(file.path('../output/SNP/CHEERS/AllPeaks/', 'Merged_condition'))
  # dir.create(file.path('../output/SNP/CHEERS/AllPeaks/', 'Individual_condition'))
  # fwrite(peak_out, file = paste0('../output/SNP/CHEERS/AllPeaks/Merged/', name, '_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
  
  for(i in names(table(meta$condition))){
    meta_sub_sub <- meta_sub[meta_sub$condition == i,]
    mat_sub_sub <- mat_sub[,rownames(meta_sub_sub)]
    mat_sub_sub_count <- rowSums(mat_sub_sub)
    peak_out <- peaks_mat_df[,c(1:3)]
    peak_out$count <- mat_sub_sub_count
    
    fwrite(peak_out, file = paste0('../output/SNP/CHEERS/AllPeaks/Merged_condition/', name, '_', i,'_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
    
    # for(j in names(table(meta$sample))){
    #   meta_sub_sub_sub <- meta_sub_sub[meta_sub_sub$sample == j,]
    #   mat_sub_sub_sub <- mat_sub_sub[,rownames(meta_sub_sub_sub)]
    #   mat_sub_sub_sub_count <- rowSums(mat_sub_sub_sub)
    #   peak_out <- peaks_mat_df[,c(1:3)]
    #   peak_out$count <- mat_sub_sub_sub_count
    #   
    #   fwrite(peak_out, file = paste0('../output/SNP/CHEERS/AllPeaks/Individual_condition/', name, '_', i,'_',j,'_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
    #   
    # }
  }
}

for(j in names(table(meta$celltype_id))){
  select_peaks_condition(j)
}
# select_peaks <- function(name){
#   peaks <- markerList[[name]]
#   peaks_gr <- peaks[,c(1,3,4,2)]
#   peaks_gr <- makeGRangesFromDataFrame(peaks_gr, keep.extra.columns = T)
#   
#   peaks_mat_sub <- subsetByOverlaps(peaks_mat, peaks_gr)
#   
#   peaks_mat_sub_df <- data.frame(chr=seqnames(peaks_mat_sub),
#                                  start=start(peaks_mat_sub),
#                                  end=end(peaks_mat_sub),
#                                  distToTSS = rowRanges(peaks_mat_sub)$distToTSS,
#                                  GC = rowRanges(peaks_mat_sub)$GC)
#   rownames(peaks_mat_sub_df) <- GRangesToString(rowRanges(peaks_mat_sub), sep = c("-", "-"))
#   
#   meta <- as.data.frame(colData(peaks_mat_sub))
#   mat <- assay(peaks_mat_sub)
#   rownames(mat) <- rownames(peaks_mat_sub_df)
#   
#   mat_sub <- mat[,rownames(meta[meta$celltype_id == name,])]
#   mat_sub_count <- rowSums(mat_sub)
#   peaks_mat_sub_df$count <- mat_sub_count
#   
#   peak_out <- peaks_mat_sub_df[,c(1:3,6)]
#   
#   dir.create(file.path('../output/SNP/CHEERS/MarkerPeaks/', name))
#   write.table(snps_out, file = paste0('../output/SNP/CHEERS/MarkerPeaks/', name, '/_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
# }
