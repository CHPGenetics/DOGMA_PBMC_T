library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
library(ggpubr)
library(cowplot)
library(patchwork)
library(DESeq2)
library(dplyr)
library(ggplot2)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code/')

df = read.xlsx('../output/colocalization/peak_snp_p2g_gene_sig_consistent_only_snp_updated.xlsx')
snp_pri <- unlist(str_split(df$SNP, '\\|'))

huang <- read.xlsx('../data/SNP/41586_2017_BFnature22969_MOESM2_ESM.xlsx', sheet = 1)
lange <- read.xlsx('../data/SNP/41588_2017_BFng3760_MOESM277_ESM.xlsx')
lange_name <- as.character(lange[1,])
lange <- lange[-c(1),]
colnames(lange) <- lange_name
huang <- huang[huang$tier2 =='No',]
lange <- lange[lange$`Confidently fine-mapped in this study?` == T,]
huang$variant <- huang$all.variant
lange$variant <- lange$`CredibleSet (minimal SNP set with >95% posterior probability)`

result <- lapply(1:nrow(huang), function(i){
    loci <- huang[i,]
    test <- intersect(unlist(str_split(loci$variant, ',')), snp_pri)
    return(ifelse(length(test) == 0, F, T))
  })
huang_loci <- huang[unlist(result),]

result <- lapply(1:nrow(lange), function(i){
  loci <- lange[i,]
  test <- intersect(unlist(str_split(loci$variant, ',')), snp_pri)
  return(ifelse(length(test) == 0, F, T))
})
lange_loci <- lange[unlist(result),]

for(i in 1:nrow(lange_loci)){
  loci <- lange_loci[i,]
  result <- unlist(lapply(str_split(df$SNP, '\\|'), function(i){
    test <- intersect(i, unlist(str_split(loci$variant, ',')))
    return(ifelse(length(test) == 0, F, T))
  }))
  genes <- unlist(str_split(df[result,]$gene_consistent, '\\|'))
  ordered_genes <- names(table(genes))[order(table(genes), decreasing = TRUE)]
  lange_loci$snps[i] <- paste0(unlist(str_split(df[result,]$SNP, '\\|')), collapse = ',')
  lange_loci$genes[i] <- paste0(ordered_genes, collapse = ',')
  }

for(i in 1:nrow(huang_loci)){
  loci <- huang_loci[i,]
  result <- unlist(lapply(str_split(df$SNP, '\\|'), function(i){
    test <- intersect(i, unlist(str_split(loci$variant, ',')))
    return(ifelse(length(test) == 0, F, T))
  }))
  genes <- unlist(str_split(df[result,]$gene_consistent, '\\|'))
  ordered_genes <- names(table(genes))[order(table(genes), decreasing = TRUE)]
  huang_loci$snps[i] <- paste0(unlist(str_split(df[result,]$SNP, '\\|')), collapse = ',')
  huang_loci$genes[i] <- paste0(ordered_genes, collapse = ',')
}

genes <- unique(c(unlist(str_split(huang_loci$genes, ',')), unlist(str_split(lange_loci$genes, ','))))

write.xlsx(huang_loci, '../data/SNP/41586_2017_BFnature22969_MOESM2_ESM_modified_sheet1.xlsx')
write.xlsx(lange_loci, '../data/SNP/41588_2017_BFng3760_MOESM277_ESM_modified.xlsx')
