library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
library(stringr)
library(ggrepel)
library(fastmatch)
set.seed(1)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/')

addArchRThreads(threads = 64)

proj <- loadArchRProject('DOGMA_filtered_multiome/')

p2g_ibd <- readRDS('../SNP/links_ibd.RDS')
p2g_ibd$value <- p2g_ibd$score
p2g <- list()
p2g[['Peak2Gene']] <- p2g_ibd

load('~/RWorkSpace/CITE-seq/Duerr/DOGMA_analysis/data/annotations.RData')

PlotAnnotation <- function(annotations, region){
  annotation <- annotations
  if (is.null(x = annotation)) {
    return(NULL)
  }
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  genes.keep <- unique(x = annotation.subset$gene_name)
  annotation.subset <- annotation[fmatch(x = annotation$gene_name, 
                                         table = genes.keep, nomatch = 0L) > 0L]
  if (length(x = annotation.subset) == 0) {
    p <- ggplot(data = data.frame())
    y_limit <- c(0, 1)
  }
  else {
    annotation_df_list <- Signac:::reformat_annotations(annotation = annotation.subset, 
                                                        start.pos = start.pos, end.pos = end.pos)
    
    peak_df <- data.frame(name = 'Gene')
    p <- ggplot(peak_df) + facet_wrap(~name, strip.position = 'right') + geom_segment(data = annotation_df_list$exons, 
                                                                                      mapping = aes_string(x = "start", y = annotation_df_list$exons$dodge, 
                                                                                                           xend = "end", yend = annotation_df_list$exons$dodge, 
                                                                                                           color = "strand"), show.legend = FALSE, size = 5) + 
      geom_segment(data = annotation_df_list$labels, mapping = aes_string(x = "start", 
                                                                          y = annotation_df_list$labels$dodge, xend = "end", 
                                                                          yend = annotation_df_list$labels$dodge, color = "strand"), 
                   show.legend = FALSE, size = 1/2)
    if (nrow(x = annotation_df_list$plus) > 0) {
      p <- p + geom_segment(data = annotation_df_list$plus, 
                            mapping = aes_string(x = "start", y = annotation_df_list$plus$dodge, 
                                                 xend = "end", yend = annotation_df_list$plus$dodge, 
                                                 color = "strand"), arrow = arrow(ends = "last", 
                                                                                  type = "open", angle = 45, length = unit(x = 0.05, 
                                                                                                                           units = "inches")), show.legend = FALSE, 
                            size = 1/2)
    }
    if (nrow(x = annotation_df_list$minus) > 0) {
      p <- p + geom_segment(data = annotation_df_list$minus, 
                            mapping = aes_string(x = "start", y = annotation_df_list$minus$dodge, 
                                                 xend = "end", yend = annotation_df_list$minus$dodge, 
                                                 color = "strand"), arrow = arrow(ends = "first", 
                                                                                  type = "open", angle = 45, length = unit(x = 0.05, 
                                                                                                                           units = "inches")), show.legend = FALSE, 
                            size = 1/2)
    }
    n_stack <- max(annotation_df_list$labels$dodge)
    annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + 
      (n_stack * 0.2)
    p <- p + geom_text(data = annotation_df_list$labels, 
                       mapping = aes_string(x = "position", y = "dodge", 
                                            label = "gene_name"), size = 3)
    y_limit <- c(0.9, n_stack + (n_stack * 0.5))
  }
  p <- p + theme_classic() + ylab("Genes") + xlab(label = paste0(chromosome, 
                                                                 " position (kb)")) + xlim(start.pos, end.pos) + ylim(y_limit) + 
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
    scale_color_manual(values = c("darkblue", "darkgreen")) + ylab("") + xlab("") + scale_x_continuous(limits = c(start(region), 
                                                                                                                  end(region)), expand = c(0, 0)) +
    theme(legend.text = element_text(size = 7)) + 
    theme_ArchR(baseSize = 7, baseLineSize = 0.4, 
                baseRectSize = 0.4) + guides(color = FALSE, 
                                             fill = FALSE) + theme(strip.text.y = element_text(size = 7, 
                                                                                               angle = 0), strip.background = element_blank()) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  return(p)
}

PlotPeak <- function(peaks, region, snp_gr){
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
  peak.intersect$group <- 'No SNP'
  peak_snp.intersect <- findOverlaps(peak.intersect, snp_gr)
  peak.intersect[queryHits(peak_snp.intersect)]$group <- 'SNP'
  
  peak.df <- as.data.frame(x = peak.intersect)
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  if (nrow(x = peak.df) > 0) {
    peak.df$start[peak.df$start < start.pos] <- start.pos
    peak.df$end[peak.df$end > end.pos] <- end.pos
    peak.df$name <- 'Peak'
    peak.plot <- ggplot(data = peak.df, aes(x = start, y = 0, xend = end, yend = 0, col = group)) + 
      facet_wrap(~name, strip.position = 'right') +
      geom_segment(size = 2, data = peak.df)+
      scale_color_manual(values = c('No SNP' = 'dodgerblue3', 'SNP' = 'firebrick'))+
      theme(legend.position = 'none')
  }
  else {
    peak.plot <- ggplot(data = peak.df)
  }
  peak.plot <- peak.plot + ylab("") + xlab("") + scale_x_continuous(limits = c(start(region), 
                                                                               end(region)), expand = c(0, 0)) +
    theme(legend.text = element_text(size = 7)) + 
    theme_ArchR(baseSize = 7, baseLineSize = 0.4, 
                baseRectSize = 0.4) + guides(color = FALSE, 
                                             fill = FALSE) + theme(strip.text.y = element_text(size = 7, 
                                                                                               angle = 0), strip.background = element_blank()) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  return(peak.plot)
}

PlotSNP <- function(peaks, region, color = "firebrick"){
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
  peak.df <- as.data.frame(x = peak.intersect)
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  if (nrow(x = peak.df) > 0) {
    peak.df$start[peak.df$start < start.pos] <- start.pos
    peak.df$end[peak.df$end > end.pos] <- end.pos
    peak.df$name <- 'SNP'
    peak.plot <- ggplot(data = peak.df, aes(x = start, y = 0, xend = end, yend = 0)) + 
      facet_wrap(~name, strip.position = 'right') +
      geom_segment(size = 2, data = peak.df, col = color)
  }
  else {
    peak.plot <- ggplot(data = peak.df)
  }
  peak.plot <- peak.plot + ylab("") + xlab("") + scale_x_continuous(limits = c(start(region), 
                                                                               end(region)), expand = c(0, 0)) +
    theme(legend.text = element_text(size = 7)) + 
    theme_ArchR(baseSize = 7, baseLineSize = 0.4, 
                baseRectSize = 0.4) + guides(color = FALSE, 
                                             fill = FALSE) + theme(strip.text.y = element_text(size = 7, 
                                                                                               angle = 0), strip.background = element_blank()) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  return(peak.plot)
}

PlotSNP_label <- function(snp_gr, peaks, region, color = "firebrick"){
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  peak.intersect <- subsetByOverlaps(x = snp_gr, ranges = region)
  peak.intersect$group <- 'No SNP'
  peak_snp.intersect <- findOverlaps(peak.intersect, peaks)
  peak.intersect[queryHits(peak_snp.intersect)]$group <- 'SNP'
  
  peak.df <- as.data.frame(x = peak.intersect)
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  if (nrow(x = peak.df) > 0) {
    peak.df$start[peak.df$start < start.pos] <- start.pos
    peak.df$end[peak.df$end > end.pos] <- end.pos
    peak.df$name <- 'SNP'
    peak.plot <- ggplot(data = peak.df, aes(x = start - 100, y = 0, xend = end + 100, yend = 0, label = SNP, col = group)) + 
      facet_wrap(~name, strip.position = 'right') +
      geom_segment(size = 2, data = peak.df)+
      geom_text_repel(angle = 90, nudge_y = 0.01, direction = 'x', segment.color = 'lightgrey')+
      scale_color_manual(values = c('No SNP' = 'dodgerblue3', 'SNP' = 'firebrick'))+
      theme(legend.position = 'none')
  }
  else {
    peak.plot <- ggplot(data = peak.df)
  }
  peak.plot <- peak.plot + ylab("") + xlab("") + scale_x_continuous(limits = c(start(region), 
                                                                               end(region)), expand = c(0, 0)) +
    theme(legend.text = element_text(size = 7)) + 
    theme_ArchR(baseSize = 7, baseLineSize = 0.4, 
                baseRectSize = 0.4) + guides(color = FALSE, 
                                             fill = FALSE) + theme(strip.text.y = element_text(size = 7, 
                                                                                               angle = 0), strip.background = element_blank()) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  return(peak.plot)
}

PlotGene_pos <- function(region, gene, color = "black"){
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  p2g_all <- Reduce(c, p2g_celltype_link)
  
  peak.intersect <- subsetByOverlaps(x = p2g_all, ranges = region, type = 'within')
  peak.df <- as.data.frame(x = peak.intersect)
  peak.df <- peak.df[!duplicated(peak.df$gene),]
  peak.df <- peak.df[peak.df$gene == gene,]
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  if (nrow(x = peak.df) > 0) {
    #peak.df$start[peak.df$start < start.pos] <- start.pos
    #peak.df$end[peak.df$end > end.pos] <- end.pos
    peak.df$name <- paste0('Position: ', gene)
    peak.plot <- ggplot(data = peak.df, aes(x = gene_pos - 100, y = 0, xend = gene_pos + 100, yend = 0)) + 
      facet_wrap(~name, strip.position = 'right') +
      geom_segment(size = 2, data = peak.df, col = color)
  }
  else {
    peak.plot <- ggplot(data = peak.df)
  }
  peak.plot <- peak.plot + ylab("") + xlab("") + scale_x_continuous(limits = c(start(region), 
                                                                               end(region)), expand = c(0, 0)) +
    theme(legend.text = element_text(size = 7)) + 
    theme_ArchR(baseSize = 7, baseLineSize = 0.4, 
                baseRectSize = 0.4) + guides(color = FALSE, 
                                             fill = FALSE) + theme(strip.text.y = element_text(size = 7, 
                                                                                               angle = 0), strip.background = element_blank()) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  return(peak.plot)
}

reform_p2g <- function(p2g){
  p2g_df <- as.data.frame(p2g)
  p2g_df$start <- p2g_df$peak_pos - 250
  p2g_df$end <- p2g_df$peak_pos + 250
  p2g_df <- p2g_df[,-c(4:5)]
  p2g_gr <- makeGRangesFromDataFrame(p2g_df, seqnames.field = 'seqnames', start.field = 'start', end.field = 'end', keep.extra.columns = T)
  return(p2g_gr)
}

PlotLink_label <- function(peaks, region, gene){
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  
  peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
  p2g_snp <- sapply(names(p2g_celltype_link), function(celltype){
    p2g <- p2g_celltype_link[[celltype]]
    p2g.intersect <- subsetByOverlaps(x = p2g, ranges = region, type = 'within')
    if (length(p2g.intersect) == 0){
      num <- 0
    }
    else{
      p2g.intersect.reform <- reform_p2g(p2g.intersect)
      p2g.snp <- subsetByOverlaps(p2g.intersect.reform, peak.intersect)
      p2g.snp <- p2g.snp[p2g.snp$gene == gene,]
      num <- length(p2g.snp)
    }
    return(num)
  })
  p2g_snp_name <- names(p2g_snp)[p2g_snp > 0]
  return(p2g_snp_name)
}

modifyTrack_label <- function(p){
  p1 <- p[[1]]$patches$plots[[4]]
  p2 <- p[[1]]$patches$plots[[1]]
  p3 <- p[[1]]$patches$plots[[2]]
  
  p1$data$y <- -p1$data$y
  p1 <- p1+coord_cartesian(ylim = c(0, 100))+
    theme(legend.position = 'none') + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  region <- p2$labels$title
  region <- gsub(' ', '', region)
  region <- gsub(':', '-', region)
  p2 <- p2+theme(strip.text.y = element_text(hjust = 0),
                 strip.background = element_blank())+ggtitle('')+labs(y = '') + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p3$labels$title <- str_split_fixed(p3$labels$title, ' : ', 2)[,2]
  p3 <- p3+theme(strip.text.y = element_blank(), 
                 strip.background = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold.italic')) + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p4 <- PlotAnnotation(annotations, region)
  p4 <- p4 + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p5 <- PlotPeak(peak_gr, region)
  p5 <- p5 +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p6 <- PlotSNP_label(snp_gr, region)
  p6 <- p6 +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p_new <- plot_spacer() + p6 + plot_spacer() + p5 + plot_spacer() + p4 + p3 + p2 + 
    plot_layout(ncol = 2, heights = c(2.1, 0.1, 1, 10), widths = c(1,10))
  return(p_new)
}

modifyTrack_label_link <- function(p){
  p1 <- p[[1]]$patches$plots[[4]]
  p2 <- p[[1]]$patches$plots[[1]]
  p3 <- p[[1]]$patches$plots[[2]]
  
  region <- p2$labels$title
  region <- gsub(' ', '', region)
  region <- gsub(':', '-', region)
  p2 <- p2+theme(strip.text.y = element_text(hjust = 0),
                 strip.background = element_blank())+ggtitle('')+labs(y = '') + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p3$labels$title <- str_split_fixed(p3$labels$title, ' : ', 2)[,2]
  p3 <- p3+theme(strip.text.y = element_blank(), 
                 strip.background = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold.italic')) + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p4 <- PlotAnnotation(annotations, region)
  p4 <- p4 + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p5 <- PlotPeak(peak_gr, region, snp_gr)
  p5 <- p5 +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p6 <- PlotSNP_label(snp_gr, peak_gr, region)
  p6 <- p6 +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p_new <- p6 + p5 + p4 + p1 +
    plot_layout(ncol = 1, heights = c(2.1, 0.1, 1, 10))
  return(p_new)
}

library(ggpubr)
library(ggsci)
modifyTrack_label_link_new <- function(p){
  p1 <- p[[1]]$patches$plots[[4]]
  p2 <- p[[1]]$patches$plots[[1]]
  p3 <- p[[1]]$patches$plots[[2]]
  
  # p1$data$y <- -p1$data$y
  # p1 <- p1+coord_cartesian(ylim = c(0, 100))+
  #   theme(legend.position = 'none') + 
  #   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  p11 <- p1 + labs(col = 'Correlation') + scale_color_gsea() + theme(legend.position = 'left')
  p12 <- as_ggplot(get_legend(p11))
  p1 <- p1 + labs(col = 'Correlation') + scale_color_gsea() + theme(legend.position = 'none', strip.text.y = element_text(hjust = 0))
  
  region <- p2$labels$title
  region <- gsub(' ', '', region)
  region <- gsub(':', '-', region)
  gene <- str_split_fixed(p3$labels$title, ' : ', 2)[,2]
  
  p2g_snp_name <- PlotLink_label(snp_gr, region, gene)
  p1$data$name[p1$data$name %in% p2g_snp_name] <- paste0(p1$data$name[p1$data$name %in% p2g_snp_name], ' *')
  
  p2 <- p2+theme(strip.text.y = element_text(hjust = 0),
                 strip.background = element_blank())+ggtitle('')+labs(y = '') + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p3$labels$title <- str_split_fixed(p3$labels$title, ' : ', 2)[,2]
  p3 <- p3+theme(strip.text.y = element_blank(), 
                 strip.background = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold.italic')) + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p4 <- PlotAnnotation(annotations, region)
  p4 <- p4 + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p5 <- PlotPeak(peak_gr, region, snp_gr)
  p5 <- p5 +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p6 <- PlotSNP_label(snp_gr, peak_gr, region)
  p6 <- p6 +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p7 <- PlotGene_pos(region, gene)
  p7 <- p7 +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  p_new <- plot_spacer() + p6 + plot_spacer() + p5 + plot_spacer() + p4 + plot_spacer() + p7 + p3 + p2 + p12 + p1 +
    plot_layout(ncol = 2, heights = c(2.6, 0.2, 0.8, 0.2, 8, 10), widths = c(1,10))
  
  return(p_new)
}

####################################
# load('DOGMA_filtered_multiome/PeakMatrix.RData')
# peak_gr = rowRanges(peaks_mat)
# save(peak_gr, file = 'DOGMA_filtered_multiome/PeakMatrix_grange.RData')
load('DOGMA_filtered_multiome/PeakMatrix_grange.RData')

snp <- read.table('../SNP/overlap_ibd_snp.bed')
colnames(snp)[4] <- 'SNP'
snp_gr <- makeGRangesFromDataFrame(snp, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', keep.extra.columns = T)

p2g_celltype_link <- list()
for(i in names(table(proj$celltype_id))){
  p2g_sub <- readRDS(paste0('DOGMA_filtered_multiome/', i, '/p2g_fdr_0.001_loop_replicate.RDS'))
  p2g_celltype_link[[i]] <- p2g_sub$Peak2GeneLinks
}

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'PTGER4', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(10, 3),
  upstream = 250000,
  downstream = 50000,
  loops = p2g,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label(p)
pdf('../../plots/visualization/Peak2Gene-Track_modify_label_PTGER4_overlap.pdf', width = 10, height = 10)
p_new
dev.off()

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'PTGER4', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(3, 10),
  upstream = 250000,
  downstream = 50000,
  loops = p2g_celltype_link,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label_link(p)
pdf('../../plots/visualization/Peak2Gene-Track_modify_label_PTGER4_overlap_link_0.001_replicate.pdf', width = 10, height = 10)
p_new
dev.off()

p_new <- modifyTrack_label_link_new(p)
pdf('../../plots/visualization/Peak2Gene-Track_modify_label_PTGER4_overlap_link_0.001_replicate_combined.pdf', width = 10, height = 10)
p_new
dev.off()

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'LNPEP', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(3, 10),
  upstream = 100000,
  downstream = 100000,
  loops = p2g_celltype_link,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label_link_new(p)
pdf('../../plots/visualization/Peak2Gene-Track_modify_label_LNPEP_overlap_link_0.001_replicate_combined.pdf', width = 10, height = 10)
p_new
dev.off()

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'IL2RA', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(3, 10),
  upstream = 60000,
  downstream = 20000,
  loops = p2g_celltype_link,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label_link_new(p)
pdf('../../plots/visualization/Peak2Gene-Track_modify_label_IL2RA_overlap_link_0.001_replicate_combined.pdf', width = 10, height = 10)
p_new
dev.off()

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'CCL20', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(3, 10),
  upstream = 20000,
  downstream = 20000,
  loops = p2g_celltype_link,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label_link_new(p)
pdf('../../plots/visualization/Peak2Gene-Track_modify_label_CCL20_overlap_link_0.001_replicate_combined.pdf', width = 10, height = 10)
p_new
dev.off()

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'CXCR5', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(3, 10),
  upstream = 20000,
  downstream = 60000,
  loops = p2g_celltype_link,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label_link_new(p)
pdf('../../plots/visualization/Peak2Gene-Track_modify_label_CXCR5_overlap_link_0.001_replicate_combined.pdf', width = 10, height = 10)
p_new
dev.off()

####################################
snp_gr_disease <- readRDS('../SNP/disease_gr.RDS')

snp_gr <- snp_gr_disease[str_detect(snp_gr_disease$Disease_merged, 'Psoriasis'),]

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'RUNX3', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(3, 10),
  upstream = 100000,
  downstream = 100000,
  loops = p2g_celltype_link,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label_link_new(p)
pdf('../../plots/visualization/Peak2Gene-Track_modify_label_RUNX3_overlap_link_0.001_replicate_combined.pdf', width = 10, height = 10)
p_new
dev.off()

snp_gr <- snp_gr_disease[str_detect(snp_gr_disease$Disease_merged, 'Behcets_disease'),]

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'KLRC4', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(3, 10),
  upstream = 20000,
  downstream = 20000,
  loops = p2g_celltype_link,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label_link_new(p)
pdf('../../plots/visualization/Peak2Gene-Track_modify_label_KLRC4_overlap_link_0.001_replicate_combined.pdf', width = 10, height = 10)
p_new
dev.off()

snp_gr <- snp_gr_disease[str_detect(snp_gr_disease$Disease_merged, 'Alopecia_areata'),]

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'IL21', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(3, 10),
  upstream = 50000,
  downstream = 50000,
  loops = p2g_celltype_link,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label_link_new(p)
pdf('../../plots/visualization/Peak2Gene-Track_modify_label_IL21_overlap_link_0.001_replicate_combined.pdf', width = 10, height = 10)
p_new
dev.off()

snp_gr <- snp_gr_disease[str_detect(snp_gr_disease$Disease_merged, 'Ankylosing_spondylitis'),]

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'GPR65', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(3, 10),
  upstream = 50000,
  downstream = 50000,
  loops = p2g_celltype_link,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label_link_new(p)
pdf('../../plots/visualization/Peak2Gene-Track_modify_label_GPR65_overlap_link_0.001_replicate_combined.pdf', width = 10, height = 10)
p_new
dev.off()

snp_gr <- snp_gr_disease[str_detect(snp_gr_disease$Disease_merged, 'Systemic_lupus_erythematosus'),]

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'CD44', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(3, 10),
  upstream = 100000,
  downstream = 50000,
  loops = p2g_celltype_link,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label_link_new(p)
pdf('../../plots/visualization/Peak2Gene-Track_modify_label_CD44_overlap_link_0.001_replicate_combined.pdf', width = 10, height = 10)
p_new
dev.off()
