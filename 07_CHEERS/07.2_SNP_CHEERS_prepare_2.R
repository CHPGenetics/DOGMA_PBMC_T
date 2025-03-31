library(data.table)
library(stringr)
library(liftOver)
library(SNPlocs.Hsapiens.dbSNP141.GRCh38)
library(rtracklayer)
library(Signac)
library(cowplot)
library(ggsci)
library(ggpubr)
library(ArchR)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')

load('../output/ArchR/DOGMA_filtered_multiome/PeakMatrix.RData')
meta <- read.csv('../output/tcell_annotated_updated.csv', row.names = 'X')

meta <- meta[meta$condition %in% c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'),]
meta$label <- factor(meta$condition, levels = c('Act_IL1B_IL23_PGE2', 'Act_IL1B_IL23'))
meta$cell_type <- meta$celltype_updated
meta$replicate <- meta$sample
cell <- meta

peaks_mat_df <- data.frame(chr=seqnames(peaks_mat),
                           start=start(peaks_mat),
                           end=end(peaks_mat),
                           distToTSS = rowRanges(peaks_mat)$distToTSS,
                           GC = rowRanges(peaks_mat)$GC)
rownames(peaks_mat_df) <- GRangesToString(rowRanges(peaks_mat), sep = c("-", "-"))

meta <- as.data.frame(colData(peaks_mat))
mat <- assay(peaks_mat)

rownames(mat) <- rownames(peaks_mat_df)
colnames(mat) <- gsub('#', '_', colnames(mat))

mat <- mat[,rownames(cell)]

rownames(meta) <- gsub('#', '_', rownames(meta))

meta <- meta[rownames(cell),]

select_peaks_condition <- function(name){
  meta_sub <- meta[meta$celltype_id == name,]
  mat_sub <- mat[,rownames(meta_sub)]
  # mat_sub_count <- rowSums(mat_sub)
  # peak_out <- peaks_mat_df[,c(1:3)]
  # peak_out$count <- mat_sub_count
  
  dir.create(file.path('../output/SNP_2/CHEERS/AllPeaks/', 'Merged_condition'))
  #dir.create(file.path('../output/SNP_2/CHEERS/AllPeaks/', 'Individual_condition'))
  # fwrite(peak_out, file = paste0('../output/SNP/CHEERS/AllPeaks/Merged/', name, '_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
  
  for(i in names(table(meta$condition))){
    meta_sub_sub <- meta_sub[meta_sub$condition == i,]
    mat_sub_sub <- mat_sub[,rownames(meta_sub_sub)]
    mat_sub_sub_count <- rowSums(mat_sub_sub)
    peak_out <- peaks_mat_df[,c(1:3)]
    peak_out$count <- mat_sub_sub_count
    
    fwrite(peak_out, file = paste0('../output/SNP_2/CHEERS/AllPeaks/Merged_condition/', name, '_', i,'_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
    
    # for(j in names(table(meta$sample))){
    #   meta_sub_sub_sub <- meta_sub_sub[meta_sub_sub$sample == j,]
    #   if(nrow(meta_sub_sub_sub) > 1){
    #     mat_sub_sub_sub <- mat_sub_sub[,rownames(meta_sub_sub_sub)]
    #     mat_sub_sub_sub_count <- rowSums(mat_sub_sub_sub)
    #     peak_out <- peaks_mat_df[,c(1:3)]
    #     peak_out$count <- mat_sub_sub_sub_count
    #     
    #     fwrite(peak_out, file = paste0('../output/SNP_2/CHEERS/AllPeaks/Individual_condition/', name, '_', i,'_',j,'_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
    #   }
    #   # else if(nrow(meta_sub_sub_sub) == 1){
    #   #   mat_sub_sub_sub <- mat_sub_sub[,rownames(meta_sub_sub_sub)]
    #   #   mat_sub_sub_sub_count <- mat_sub_sub_sub
    #   #   peak_out <- peaks_mat_df[,c(1:3)]
    #   #   peak_out$count <- mat_sub_sub_sub_count
    #   #   
    #   #   fwrite(peak_out, file = paste0('../output/SNP_2/CHEERS/AllPeaks/Individual_condition/', name, '_', i,'_',j,'_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
    #   # }
    # }
  }
}

for(j in names(table(meta$celltype_id))){
  select_peaks_condition(j)
}

set.seed(2023)
sample_select <- sample(names(table(meta$sample)), 12)
# [1] "SB775567" "SB775496" "SB775483" "SB774131" "SB774120" "SB775522" "SB775449" "SB773931"
# [9] "SB775505" "SB770848" "SB775465" "SB775598"
meta <- meta[meta$sample %in% sample_select,]

select_peaks_condition_downsample <- function(name){
  meta_sub <- meta[meta$celltype_id == name,]
  mat_sub <- mat[,rownames(meta_sub)]
  # mat_sub_count <- rowSums(mat_sub)
  # peak_out <- peaks_mat_df[,c(1:3)]
  # peak_out$count <- mat_sub_count
  
  #dir.create(file.path('../output/SNP_2/CHEERS/AllPeaks/', 'Merged_condition'))
  dir.create(file.path('../output/SNP_2/CHEERS/AllPeaks/', 'Individual_condition'))
  # fwrite(peak_out, file = paste0('../output/SNP/CHEERS/AllPeaks/Merged/', name, '_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
  
  for(i in names(table(meta$condition))){
    meta_sub_sub <- meta_sub[meta_sub$condition == i,]
    mat_sub_sub <- mat_sub[,rownames(meta_sub_sub)]
    mat_sub_sub_count <- rowSums(mat_sub_sub)
    peak_out <- peaks_mat_df[,c(1:3)]
    peak_out$count <- mat_sub_sub_count
    
    #fwrite(peak_out, file = paste0('../output/SNP_2/CHEERS/AllPeaks/Merged_condition/', name, '_', i,'_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
    
    for(j in names(table(meta$sample))){
      meta_sub_sub_sub <- meta_sub_sub[meta_sub_sub$sample == j,]
      if(nrow(meta_sub_sub_sub) > 1){
        mat_sub_sub_sub <- mat_sub_sub[,rownames(meta_sub_sub_sub)]
        mat_sub_sub_sub_count <- rowSums(mat_sub_sub_sub)
        peak_out <- peaks_mat_df[,c(1:3)]
        peak_out$count <- mat_sub_sub_sub_count
        
        fwrite(peak_out, file = paste0('../output/SNP_2/CHEERS/AllPeaks/Individual_condition/', name, '_', i,'_',j,'_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
      }
      # else if(nrow(meta_sub_sub_sub) == 1){
      #   mat_sub_sub_sub <- mat_sub_sub[,rownames(meta_sub_sub_sub)]
      #   mat_sub_sub_sub_count <- mat_sub_sub_sub
      #   peak_out <- peaks_mat_df[,c(1:3)]
      #   peak_out$count <- mat_sub_sub_sub_count
      #   
      #   fwrite(peak_out, file = paste0('../output/SNP_2/CHEERS/AllPeaks/Individual_condition/', name, '_', i,'_',j,'_ReadsInPeaks.txt'), quote = F, row.names = F, col.names = T, sep = "\t")
      # }
    }
  }
}

for(j in names(table(meta$celltype_id))){
  select_peaks_condition_downsample(j)
}
