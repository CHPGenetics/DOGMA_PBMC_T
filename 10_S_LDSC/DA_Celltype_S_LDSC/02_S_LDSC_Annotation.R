rm(list=ls())
gc()

# Load Packages
library(data.table)
library(readxl)
library(tidyverse)
library(liftOver)
library(rtracklayer)

# Parameters
REF_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Input/"
WORK_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Input/Pseudo_Celltype_DA/"
LIFT_PATH <- "/ix1/wchen/Shiyue/References/LiftOver/"

# Load Data
all_clusters <- readRDS(paste0(REF_PATH, "DOGMA_Celltypes.rds"))

# Peak3 Ref
DA_ref <- data.frame()
for(celltype in all_clusters){
  
  DA_ref_i <- fread(paste0(WORK_PATH, "DA_", gsub(" ", "_", celltype), ".txt"))
  DA_ref <- rbind(DA_ref, DA_ref_i)
}

# hg38 to hg19
Peak <- data.frame(Peak = unique(DA_ref$Peak))
Peak <- Peak %>%
  separate(Peak, into = c("seqnames", "start", "end"), sep = "-", remove = FALSE) %>%
  mutate(start = as.numeric(start)) %>%
  mutate(end = as.numeric(end)) %>%
  filter(seqnames %in% paste0("chr", 1:22))

## Liftover
cur <- makeGRangesFromDataFrame(Peak, keep.extra.columns = T)
ch <- import.chain(paste0(LIFT_PATH, 'hg38ToHg19.over.chain'))
Peaks <- GRanges(
  seqnames = Rle(Peak$seqnames),
  ranges = IRanges(Peak$start, Peak$end),
  mcol = Peak$Peak
)

lifted <- liftOver(Peaks, ch)
flattened_lifted <- unlist(lifted, use.names = TRUE)
lifted_df <- as.data.frame(flattened_lifted)
colnames(lifted_df)[which(colnames(lifted_df) == "mcol")] <- "original_peak_name"
lifted_df <- lifted_df %>%
  group_by(original_peak_name) %>%
  filter(width == max(width)) %>%
  ungroup()

lifted_df <- lifted_df %>%
  distinct(original_peak_name, .keep_all = TRUE) %>%
  mutate(seqnames = gsub("chr", "", seqnames)) %>%
  mutate(seqnames = as.numeric(seqnames)) %>%
  mutate(start = as.numeric(start)) %>%
  mutate(end = as.numeric(end))

Peak3 <- data.frame(GENE = lifted_df$original_peak_name,
                    CHR = lifted_df$seqnames,
                    START = lifted_df$start,
                    END = lifted_df$end) %>%
  arrange(CHR, START, END)
write.table(Peak3, paste0(REF_PATH, "Peak3_coord_all.txt"),
            sep = "\t", quote = F, col.names = T, row.names = F)
write.table(Peak3$GENE, paste0(REF_PATH, "Peak3_list.txt"),
            sep = "\t", quote = F, col.names = F, row.names = F)

n_peak <- round(nrow(Peak3)/5)
for(celltype in all_clusters){
  
  DA <- fread(paste0(WORK_PATH, "DA_", gsub(" ", "_", celltype), ".txt")) %>%
    filter(stringr::str_extract(Peak, "^[^-]+") %in% paste0("chr", 1:22)) %>%
    filter(Peak %in% Peak3$GENE)
  DA <- DA %>%
    arrange(p_val)
  
  DA_sub <- DA[1:n_peak,]
  write.table(DA_sub$Peak,
              file = paste0(WORK_PATH, "Peak_list_", gsub(" ", "_", celltype), ".txt"),
              sep = "\t", quote = F, col.names = F, row.names = F)
}

# ldcts
LDSC_PATH <- "/ix1/wchen/Shiyue/Projects/2023_07_DOGMA_Revision/Results/02_S_LDSC/Output/"
cluster_file <- gsub(" ", "_", all_clusters)
path_list <- paste0(LDSC_PATH, "Pseudo_Celltype_DA/", cluster_file, "_", ",", LDSC_PATH, "Annotation_Peak3/ref_")
ldcts <- data.frame(cluster_file, path_list)
write.table(ldcts, row.names = F, col.names = F, quote = F, sep = "\t",
            file = paste0(LDSC_PATH, "Pseudo_Celltype_DA/Result/ldcts"))