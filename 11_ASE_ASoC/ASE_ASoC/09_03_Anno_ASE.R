rm(list=ls())
gc()

# Load Packages
library(bigreadr)
library(dplyr)
library(stringr)
library(ggplot2)
library(glue)

### 0. format for anovar
ASoC_count <- fread2("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASoC_Count/ASoC_counts.txt")
ASE_count <- fread2("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASE_Count/ASE_counts.txt")

anovar_asoc <- ASoC_count[,c(1, 2, 2, 4, 5, 3)]
anovar_asoc$contig <- gsub("chr", "", anovar_asoc$contig)
anovar_asoc <- anovar_asoc[!duplicated(anovar_asoc$variantID),]
anovar_ase <- ASE_count[,c(1, 2, 2, 4, 5, 3)]
anovar_ase$contig <- gsub("chr", "", anovar_ase$contig)
anovar_ase <- anovar_ase[!duplicated(anovar_ase$variantID),]

#
fwrite2(anovar_asoc, file = "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/anovar_asoc.avinput",
        sep = "\t", col.names = F)
fwrite2(anovar_ase, file = "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/anovar_ase.avinput",
        sep = "\t", col.names = F)
# #
# anovar_asoc_refGene <- fread2("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/anovar_asoc_refGene.variant_function")
# anovar_ase_refGene <- fread2("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/anovar_ase_refGene.variant_function")

### 1. extract gene list for variants
ASE_gene_df <- str_split(anovar_ase_refGene$V2, "\\),", simplify = T) %>% as.data.frame()
ASE_gene_df2 <- lapply(1:ncol(ASE_gene_df), function(x){
  str_split(ASE_gene_df[,x], "\\(", simplify = T) %>% as.data.frame()
}) %>% Reduce("cbind", .)
ASE_gene_df3 <- lapply(1:ncol(ASE_gene_df2), function(x){
  
  xx <- ASE_gene_df2[,x]
  xx[grep("\\)", xx)] <- NA
  xx[grep("\\=", xx)] <- NA
  xx[grep("\\:", xx)] <- NA
  xx[xx == ""] <- NA
  return(xx)  
  
}) %>% Reduce("cbind", .)
ASE_gene_df4 <- lapply(1:ncol(ASE_gene_df3), function(x){
  
  str_split(ASE_gene_df3[,x], ",", simplify = T) %>% as.data.frame()
  
}) %>% Reduce("cbind", .)
ASE_gene_list <- lapply(1:nrow(ASE_gene_df4), function(x){
  
  xx <- ASE_gene_df4[x,]
  xx <- xx[!is.na(xx)]
  xx <- xx[xx != ""]
  
})

##
ASOC_gene_df <- str_split(anovar_asoc_refGene$V2, "\\),", simplify = T) %>% as.data.frame()
ASOC_gene_df2 <- lapply(1:ncol(ASOC_gene_df), function(x){
  str_split(ASOC_gene_df[,x], "\\(", simplify = T) %>% as.data.frame()
}) %>% Reduce("cbind", .)
ASOC_gene_df3 <- lapply(1:ncol(ASOC_gene_df2), function(x){
  
  xx <- ASOC_gene_df2[,x]
  xx[grep("\\)", xx)] <- NA
  xx[grep("\\=", xx)] <- NA
  xx[grep("\\:", xx)] <- NA
  xx[xx == ""] <- NA
  return(xx)  
  
}) %>% Reduce("cbind", .)
ASOC_gene_df4 <- lapply(1:ncol(ASOC_gene_df3), function(x){
  
  str_split(ASOC_gene_df3[,x], ",", simplify = T) %>% as.data.frame()
  
}) %>% Reduce("cbind", .)
ASOC_gene_list <- lapply(1:nrow(ASOC_gene_df4), function(x){
  
  xx <- ASOC_gene_df4[x,]
  xx <- xx[!is.na(xx)]
  xx <- xx[xx != ""]
  
})

### 2. get gene infor from ENSEMBL
all_gene <- c(unlist(ASOC_gene_list), unlist(ASE_gene_list)) %>% unique
all_gene_ens <- all_gene[grep("^ENSG", all_gene)]
all_gene_symbol <- setdiff(all_gene, all_gene_ens)
##
library(biomaRt)
library(httr)

set_config(config(ssl_verifypeer = 0L))
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl",
                      GRCh = 38)

gene_bp1 <- getBM(attributes = c("chromosome_name", 
                                 "transcription_start_site",
                                 "external_gene_name",
                                 "ensembl_gene_id"), 
                  filters = "external_gene_name", 
                  values = all_gene_symbol, 
                  mart = ensembl)
gene_bp2 <- getBM(attributes = c("chromosome_name", 
                                 "transcription_start_site",
                                 "external_gene_name",
                                 "ensembl_gene_id"), 
                  filters = "ensembl_gene_id", 
                  values = all_gene_ens, 
                  mart = ensembl)

gene_bp <- rbind(gene_bp1, gene_bp2)
gene_bpf <- subset(gene_bp, gene_bp$chromosome_name %in% 1:22)
gene_bpf$chromosome_name <- as.integer(gene_bpf$chromosome_name)

gene_bpf_ens <- gene_bpf[match(all_gene_ens, gene_bpf$ensembl_gene_id),]
gene_bpf_ens <- gene_bpf_ens[!is.na(gene_bpf_ens$ensembl_gene_id),]
gene_bpf_ens <- gene_bpf_ens[!duplicated(gene_bpf_ens$ensembl_gene_id), ]

gene_bpf_symbol <- gene_bpf[match(all_gene_symbol, gene_bpf$external_gene_name),]
gene_bpf_symbol <- gene_bpf_symbol[!is.na(gene_bpf_symbol$external_gene_name),]
gene_bpf_symbol <- gene_bpf_symbol[!duplicated(gene_bpf_symbol$external_gene_name), ]

colnames(gene_bpf_ens) <- colnames(gene_bpf_symbol) <- 
  c("chr", "tss", "gene1", "gene2")
rownames(gene_bpf_ens) <- rownames(gene_bpf_symbol) <- NULL

saveRDS(gene_bpf_ens, file = "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/gene_bpf_ens.rds")
saveRDS(gene_bpf_symbol, file = "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/gene_bpf_symbol.rds")

## format gene infor
ASOC_genebp_df <- lapply(1:length(ASOC_gene_list), function(x){
  
  gene_listx <- ASOC_gene_list[[x]]
  anno_dfx <- rbind(gene_bpf_ens[match(intersect(all_gene_ens, gene_listx),
                                       gene_bpf_ens$gene2),],
                    gene_bpf_symbol[match(intersect(all_gene_symbol, gene_listx),
                                          gene_bpf_symbol$gene1),])
  anno_dfx$variantID <- anovar_asoc_refGene$V8[x]
  anno_dfx$Pos <- anovar_asoc_refGene$V4[x]
  return(anno_dfx)
}) %>% Reduce("rbind", .)


ASE_genebp_df <- lapply(1:length(ASE_gene_list), function(x){
  
  gene_listx <- ASE_gene_list[[x]]
  anno_dfx <- rbind(gene_bpf_ens[match(intersect(all_gene_ens, gene_listx),
                                       gene_bpf_ens$gene2),],
                    gene_bpf_symbol[match(intersect(all_gene_symbol, gene_listx),
                                          gene_bpf_symbol$gene1),])
  anno_dfx$variantID <- anovar_ase_refGene$V8[x]
  anno_dfx$Pos <- anovar_ase_refGene$V4[x]
  return(anno_dfx)
}) %>% Reduce("rbind", .)

ASOC_genebp_df$dis2tss <- ASOC_genebp_df$Pos - ASOC_genebp_df$tss
ASOC_genebp_dff <- ASOC_genebp_df[order(abs(ASOC_genebp_df$dis2tss), decreasing = F),]
ASOC_genebp_dff <- ASOC_genebp_dff[!duplicated(ASOC_genebp_dff$variantID),]
ASE_genebp_df$dis2tss <- ASE_genebp_df$Pos - ASE_genebp_df$tss
ASE_genebp_dff <- ASE_genebp_df[order(abs(ASE_genebp_df$dis2tss), decreasing = F),]
ASE_genebp_dff <- ASE_genebp_dff[!duplicated(ASE_genebp_dff$variantID),]

saveRDS(ASOC_genebp_df, file = "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/ASOC_genebp_df.rds")
saveRDS(ASE_genebp_df, file = "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/ASE_genebp_df.rds")

anovar_asoc_refGene$gene_nearest <- ASOC_genebp_dff$gene1[match(anovar_asoc_refGene$V8, 
                                                                ASOC_genebp_dff$variantID)]
anovar_asoc_refGene$tss <- ASOC_genebp_dff$tss[match(anovar_asoc_refGene$V8, 
                                                     ASOC_genebp_dff$variantID)]
##
anovar_ase_refGene$gene_nearest <- ASE_genebp_dff$gene1[match(anovar_ase_refGene$V8, 
                                                              ASE_genebp_dff$variantID)]
anovar_ase_refGene$tss <- ASE_genebp_dff$tss[match(anovar_ase_refGene$V8, 
                                                   ASE_genebp_dff$variantID)]
saveRDS(anovar_ase_refGene, file = "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/anovar_ase_refGene.rds")
saveRDS(anovar_ase_refGene, file = "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/anovar_ase_refGene.rds")

### 3. add gene infor to count data
ASoC_count <- fread2("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASoC_Count/ASoC_counts.txt")
ASE_count <- fread2("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/04_ASE_Count/ASE_counts.txt")

ASoC_count$ref_anno <- anovar_asoc_refGene$V1[match(ASoC_count$variantID, anovar_asoc_refGene$V8)]
ASoC_count$ref_gene <- anovar_asoc_refGene$V2[match(ASoC_count$variantID, anovar_asoc_refGene$V8)]
ASoC_count$gene_nearest <- anovar_asoc_refGene$gene_nearest[match(ASoC_count$variantID, anovar_asoc_refGene$V8)]
ASoC_count$tss <- anovar_asoc_refGene$tss[match(ASoC_count$variantID, anovar_asoc_refGene$V8)]
ASoC_count$dist_to_tss <- ASoC_count$position - ASoC_count$tss

ASE_count$ref_anno <- anovar_ase_refGene$V1[match(ASE_count$variantID, anovar_ase_refGene$V8)]
ASE_count$ref_gene <- anovar_ase_refGene$V2[match(ASE_count$variantID, anovar_ase_refGene$V8)]
ASE_count$gene_nearest <- anovar_ase_refGene$gene_nearest[match(ASE_count$variantID, anovar_ase_refGene$V8)]
ASE_count$tss <- anovar_ase_refGene$tss[match(ASE_count$variantID, anovar_ase_refGene$V8)]
ASE_count$dist_to_tss <- ASE_count$position - ASE_count$tss

saveRDS(ASoC_count, file = "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/ASoC_count.rds")
saveRDS(ASE_count, file = "/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/ASE_count.rds")

### 4. plot
library(bigreadr)
library(dplyr)
library(stringr)
library(ggplot2)
library(glue)
library(cols4all)
library(patchwork)
#
group_use_anno <- c4a("Set1", 4)
col_use_anno <- c4a("rainbow", 12)
target_gene <- c("IL23R", "IL1R1", "PTGER4")
##
ASoC_count <- readRDS("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/ASoC_count.rds")
ASE_count <- readRDS("/ix1/wchen/Shiyue/Projects/2024_05_Allelic_Imbalance_Method/ASE/Result/05_Plots/Anno/ASE_count.rds")
ASoC_countf <- subset(ASoC_count, ASoC_count$p_value_adj < 0.05)
ASE_countf <- subset(ASE_count, ASE_count$p_value_adj < 0.05)

## density plot for distance
dens_asoc <- ggplot(ASoC_countf) + 
  geom_density(aes(x = dist_to_tss, 
                   color = condition, 
                   fill = condition), alpha = 0.3) + 
  facet_wrap(~ASoC_countf$celltype, scales = "free_y") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11, color = "black"))
#
dens_ase <- ggplot(ASE_countf) + 
  geom_density(aes(x = dist_to_tss, 
                   color = condition, 
                   fill = condition), alpha = 0.3) + 
  facet_wrap(~ASE_countf$celltype, scales = "free_y") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11, color = "black"))
## box plot for distance
box_plt_asoc <- ggplot(ASoC_countf) + 
  geom_boxplot(aes(x = condition, 
                   y = dist_to_tss,
                   color = condition), 
               outliers = F, width = 0.6, size = 1) + 
  xlab("") + ylab("Distance to nearest gene") +
  ggtitle("ASoC") +
  scale_color_manual(values = group_use_anno) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        title = element_text(size = 15, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none")
#
box_plt_ase <- ggplot(ASE_countf) + 
  geom_boxplot(aes(x = condition, 
                   y = dist_to_tss,
                   color = condition), 
               outliers = F, width = 0.6, size = 1) + 
  xlab("") + ylab("Distance to nearest gene") +
  ggtitle("ASE") +
  scale_color_manual(values = group_use_anno) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        title = element_text(size = 15, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11, color = "black"))
#
box_plt_asocg <- ggplot(ASoC_countf) + 
  geom_boxplot(aes(x = condition, 
                   y = dist_to_tss,
                   color = condition), 
               outliers = F, width = 0.6, size = 1) + 
  xlab("") + ylab("Distance to nearest gene") +
  ggtitle("ASoC") +
  facet_wrap(~ASoC_countf$celltype, scales = "free_y") + 
  scale_color_manual(values = group_use_anno) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        title = element_text(size = 15, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11, color = "black"))
#
box_plt_aseg <- ggplot(ASE_countf) + 
  geom_boxplot(aes(x = condition, 
                   y = dist_to_tss,
                   color = condition), 
               outliers = F, width = 0.6, size = 1) + 
  xlab("") + ylab("Distance to nearest gene") +
  ggtitle("ASE") +
  facet_wrap(~ASE_countf$celltype, scales = "free_y") + 
  scale_color_manual(values = group_use_anno) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        title = element_text(size = 15, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11, color = "black"))

## box plot for annotation
ASoC_countf$ref_anno[grep(";", ASoC_countf$ref_anno)] <- "Mixed"
ptab_asco_group <- table(ASoC_countf$ref_anno, 
                         ASoC_countf$condition) %>% 
  prop.table(.,2) %>% as.data.frame()
colnames(ptab_asco_group) <- c("Relation_to_near_genes", "Group", "Proportion")
ptab_asco_group$Proportion <- ptab_asco_group$Proportion * 100
ptab_asco_group$Relation_to_near_genes <- factor(ptab_asco_group$Relation_to_near_genes,
                                                 levels = c("exonic",
                                                            "splicing",
                                                            "ncRNA_exonic",
                                                            "ncRNA_splicing",
                                                            "UTR3", "UTR5",
                                                            "intronic",
                                                            "ncRNA_intronic",
                                                            "upstream", "downstream",
                                                            "intergenic", "Mixed"))
#
bar_plt_asoc <- ggplot(ptab_asco_group) + 
  geom_bar(aes(x = Group, y = Proportion, fill = Relation_to_near_genes),
           position = "stack", stat="identity", width = 0.8) + 
  scale_fill_manual(values = col_use_anno) +
  xlab("") + ylab("Proportion(%)") + 
  ggtitle("ASoC") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        title = element_text(size = 15, face = "bold"),
        legend.position = "none",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11, color = "black"))

#
ASE_countf$ref_anno[grep(";", ASE_countf$ref_anno)] <- "Mixed"
ptab_ase_group <- table(ASE_countf$ref_anno, 
                        ASE_countf$condition) %>% 
  prop.table(.,2) %>% as.data.frame()
colnames(ptab_ase_group) <- c("Relation_to_near_genes", "Group", "Proportion")
ptab_ase_group$Proportion <- ptab_ase_group$Proportion * 100
ptab_ase_group$Relation_to_near_genes <- factor(ptab_ase_group$Relation_to_near_genes,
                                                levels = c("exonic",
                                                           "splicing",
                                                           "ncRNA_exonic",
                                                           "ncRNA_splicing",
                                                           "UTR3", "UTR5",
                                                           "intronic",
                                                           "ncRNA_intronic",
                                                           "upstream", "downstream",
                                                           "intergenic", "Mixed"))
#
bar_plt_ase <- ggplot(ptab_ase_group) + 
  geom_bar(aes(x = Group, y = Proportion, fill = Relation_to_near_genes),
           position = "stack", stat="identity", width = 0.8) + 
  scale_fill_manual(values = col_use_anno) +
  xlab("") + ylab("Proportion(%)") + 
  ggtitle("ASE") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11, color = "black"))
