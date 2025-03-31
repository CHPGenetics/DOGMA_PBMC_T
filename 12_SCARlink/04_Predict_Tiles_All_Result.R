rm(list=ls())
gc()

# Load Packages
library(tidyverse)
library(Seurat)
library(data.table)
library(ggplot2)
library(openxlsx)


result <- fread("/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Result/All/gene_linked_tiles_celltype_updated.csv.gz")

result_f <- result %>% filter(FDR < 0.05)
write.xlsx(result_f, file = "/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Result/All/gene_linked_tiles.xlsx")

summary <- data.frame()
for(i in unique(result_f$celltype_updated)){
  
  result_s <- result_f %>% filter(celltype_updated == i)
  info <- c(celltype = i, tile_num = nrow(result_s), gene_num = length(unique(result_s$gene)))
  summary <- rbind(summary, info)
}

# Validate P2G
levels = c('CD4+ Naive (Resting)', 'CD4+ Naive (Activated)',
           'CD4+ Regulatory (Resting)', 'CD4+ Regulatory (Activated)',
           'CD4+ Memory (Resting) - Th1', 'CD4+ Memory (Activated) - Th1',
           'CD4+ Memory (Resting) - Th17', 'CD4+ Memory (Activated) - Th17',
           'CD4+ Memory (Resting) - Tfh', 'CD4+ Memory (Activated) - Tfh',
           'CD4+ Memory (Resting) - Other', 'CD4+ Memory (Activated) - Other',
           'CD8+ Naive (Resting)', 'CD8+ Naive (Activated)',
           'CD8+ Regulatory',
           'CD8+ Memory (Resting)', 'CD8+ Memory (Activated)',
           'MAITs (Resting)', 'MAITs (Activated)', 'Gamma Delta')

summary <- data.frame()
colocal <- data.frame()
for(i in unique(result_f$celltype_updated)){

  position <- match(i, levels)
  formatted_position <- sprintf("%02d", position)
  p2g <- readRDS(paste0('/ix1/wchen/zhongli/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/DOGMA_filtered_multiome/', 
                        formatted_position, "_", i, '/p2g_fdr_0.05_loop_replicate_peak_gene.RDS'))
  
  p2g <- data.frame(p2g$Peak2GeneLinks)
  p2g <- p2g %>% filter(FDR < 0.05)
  result_s <- result_f %>% filter(celltype_updated == i)
  
  genes <- intersect(result_s$gene, p2g$gene)
  result_s <- result_s %>% filter(gene %in% genes)
  p2g <- p2g %>% filter(gene %in% genes)
  
  colocal_s <- p2g %>%
    rowwise() %>%
    do({
      peak <- .
      matching_tiles <- result_s %>%
        filter(
          chr == peak$seqnames,
          (start >= peak$start & end <= peak$end) |  # tile 完全包含在 peak 中
            (start <= peak$start & end >= peak$start) |  # tile 穿过 peak 的 start
            (start <= peak$end & end >= peak$end)  # tile 穿过 peak 的 end
        )
      if (nrow(matching_tiles) > 0) {
        cbind(matching_tiles, peak_start = peak$start, peak_end = peak$end)
      } else {
        data.frame()
      }
    }) %>%
    bind_rows()
  
  colocal_ss <- colocal_s %>% select(peak_start, peak_end, gene) %>% distinct()
  
  info <- c(celltype = i, num = nrow(colocal_ss))
  summary <- rbind(summary, info)
  colocal <- rbind(colocal, colocal_s)
}
saveRDS(colocal, "/ix1/rduerr/shared/rduerr_wchen/Shiyue/Chen_Files/2023_07_DOGMA_Revision/SCARlink/Result/All/coloc_tile_p2g.rds")

summary[,2] <- as.numeric(summary[,2])
mean(summary[,2])

# Check PTGER4
colocal <- colocal %>% ungroup()
colocal_s <- colocal %>% filter(gene == "CCL20")

a <- result %>% select(gene, `Spearman corr`)
a <- a %>% distinct()


