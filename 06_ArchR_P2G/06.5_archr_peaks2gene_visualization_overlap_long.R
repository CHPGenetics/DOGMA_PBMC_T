library(ArchR)
library(stringr)
library(Seurat)
library(Signac)
library(openxlsx)
library(stringr)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
library(EnsDb.Hsapiens.v86)
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

# p2g <- getPeak2GeneLinks(
#   ArchRProj = proj,
#   corCutOff = 0,
#   FDRCutOff= 0.05,
#   resolution = 1
# )

# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "hg38"

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

PlotPeak <- function(peaks, region, color = "dodgerblue3"){
  # assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  # if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
  #   stop("The requested assay is not a ChromatinAssay.")
  # }
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  # if (is.null(x = peaks)) {
  #   peaks <- granges(x = object[[assay]])
  #   md <- object[[assay]][[]]
  #   mcols(x = peaks) <- md
  # }
  peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
  peak.df <- as.data.frame(x = peak.intersect)
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  if (nrow(x = peak.df) > 0) {
    # if (!is.null(x = group.by)) {
    #   if (!(group.by %in% colnames(x = peak.df))) {
    #     warning("Requested grouping variable not found")
    #     group.by <- NULL
    #   }
    # }
    peak.df$start[peak.df$start < start.pos] <- start.pos
    peak.df$end[peak.df$end > end.pos] <- end.pos
    peak.df$name <- 'Peak'
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
  # if (is.null(x = group.by)) {
  #   peak.plot <- peak.plot + scale_color_manual(values = color) + 
  #     theme(legend.position = "none")
  # }
  return(peak.plot)
}

PlotSNP <- function(peaks, region, color = "firebrick"){
  # assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  # if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
  #   stop("The requested assay is not a ChromatinAssay.")
  # }
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  # if (is.null(x = peaks)) {
  #   peaks <- granges(x = object[[assay]])
  #   md <- object[[assay]][[]]
  #   mcols(x = peaks) <- md
  # }
  peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
  peak.df <- as.data.frame(x = peak.intersect)
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  if (nrow(x = peak.df) > 0) {
    # if (!is.null(x = group.by)) {
    #   if (!(group.by %in% colnames(x = peak.df))) {
    #     warning("Requested grouping variable not found")
    #     group.by <- NULL
    #   }
    # }
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
  # if (is.null(x = group.by)) {
  #   peak.plot <- peak.plot + scale_color_manual(values = color) + 
  #     theme(legend.position = "none")
  # }
  return(peak.plot)
}

PlotSNP_label <- function(peaks, region, color = "firebrick"){
  # assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  # if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
  #   stop("The requested assay is not a ChromatinAssay.")
  # }
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  # if (is.null(x = peaks)) {
  #   peaks <- granges(x = object[[assay]])
  #   md <- object[[assay]][[]]
  #   mcols(x = peaks) <- md
  # }
  peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
  peak.df <- as.data.frame(x = peak.intersect)
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  if (nrow(x = peak.df) > 0) {
    # if (!is.null(x = group.by)) {
    #   if (!(group.by %in% colnames(x = peak.df))) {
    #     warning("Requested grouping variable not found")
    #     group.by <- NULL
    #   }
    # }
    peak.df$start[peak.df$start < start.pos] <- start.pos
    peak.df$end[peak.df$end > end.pos] <- end.pos
    peak.df$name <- 'SNP'
    peak.plot <- ggplot(data = peak.df, aes(x = start, y = 0, xend = end, yend = 0, label = SNP)) + 
      facet_wrap(~name, strip.position = 'right') +
      geom_segment(size = 2, data = peak.df, col = color)+
      geom_text_repel(angle = 90, nudge_y = 0.01, direction = 'x', segment.color = 'gray')
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
  # if (is.null(x = group.by)) {
  #   peak.plot <- peak.plot + scale_color_manual(values = color) + 
  #     theme(legend.position = "none")
  # }
  return(peak.plot)
}

peak <- read.table('../SNP/overlap_ibd_peak.bed')
peak$chr_pos <- paste(peak$V1, peak$V2, peak$V3, sep = '-')
peak <- peak[!duplicated(peak$chr_pos),]
peak_gr <- makeGRangesFromDataFrame(peak, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')

snp <- read.table('../SNP/overlap_ibd_snp.bed')
colnames(snp)[4] <- 'SNP'
snp_gr <- makeGRangesFromDataFrame(snp, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', keep.extra.columns = T)

####################################
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'PTGER4', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(10, 3),
  upstream = 350000,
  downstream = 50000,
  loops = p2g,
  useMatrix = 'GeneExpressionMatrix'
)

# plotPDF(p, name = "Peak2Gene-Track", width = 12, height = 12, ArchRProj = proj, addDOC = FALSE)

# modifyTrack <- function(p){
#   p1 <- p[[1]]$patches$plots[[4]]
#   p2 <- p[[1]]$patches$plots[[1]]
#   p3 <- p[[1]]$patches$plots[[2]]
#   
#   p1$data$y <- -p1$data$y
#   p1 <- p1+coord_cartesian(ylim = c(0, 100))+
#     theme(legend.position = 'none') + 
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#   
#   region <- p2$labels$title
#   region <- gsub(' ', '', region)
#   region <- gsub(':', '-', region)
#   p2 <- p2+theme(strip.text.y = element_text(hjust = 0),
#                  strip.background = element_blank())+ggtitle('')+labs(y = '') + 
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#   
#   p3$labels$title <- str_split_fixed(p3$labels$title, ' : ', 2)[,2]
#   p3 <- p3+theme(strip.text.y = element_blank(), 
#                  strip.background = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold.italic')) + 
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#   
#   p4 <- PlotAnnotation(annotations, region)
#   p4 <- p4 + 
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#   
#   p5 <- PlotPeak(peak_gr, region)
#   p5 <- p5 +
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#   
#   p6 <- PlotSNP(snp_gr, region)
#   p6 <- p6 +
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#   
#   p_new <- plot_spacer() + p1 + plot_spacer() + p6 + plot_spacer() + p5 + plot_spacer() + p4 + p3 + p2 + 
#     plot_layout(ncol = 2, heights = c(2, 0.1, 0.1, 1, 10), widths = c(1,10))
#   return(p_new)
# }
# 
# p_new <- modifyTrack(p)
# pdf('DOGMA_filtered_multiome/Plots/Peak2Gene-Track_modify.pdf', width = 12, height = 12)
# p_new
# dev.off()

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

p_new <- modifyTrack_label(p)
pdf('DOGMA_filtered_multiome/Plots/Peak2Gene-Track_modify_label_PTGER4_overlap_long.pdf', width = 10, height = 10)
p_new
dev.off()

modifyTrack_label_link <- function(p){
  p1 <- p[[1]]$patches$plots[[4]]
  p2 <- p[[1]]$patches$plots[[1]]
  p3 <- p[[1]]$patches$plots[[2]]
  
  # p1$data$y <- -p1$data$y
  # p1 <- p1+coord_cartesian(ylim = c(0, 100))+
  #   theme(legend.position = 'none') + 
  #   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
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
  
  p_new <- p6 + p5 + p4 + p1 +
    plot_layout(ncol = 1, heights = c(2.1, 0.1, 1, 10))
  return(p_new)
}

p2g_celltype_link <- list()
for(i in names(table(proj$celltype_id))){
  p2g_sub <- readRDS(paste0('DOGMA_filtered_multiome/', i, '/p2g_fdr_0.05_loop.RDS'))
  p2g_celltype_link[[i]] <- p2g_sub$Peak2GeneLinks
}

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'PTGER4', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(3, 10),
  upstream = 350000,
  downstream = 50000,
  loops = p2g_celltype_link,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label_link(p)
pdf('DOGMA_filtered_multiome/Plots/Peak2Gene-Track_modify_label_PTGER4_overlap_link_long.pdf', width = 10, height = 10)
p_new
dev.off()

############
p2g_celltype_link <- list()
for(i in names(table(proj$celltype_id))){
  p2g_sub <- readRDS(paste0('DOGMA_filtered_multiome/', i, '/p2g_fdr_0.01_loop.RDS'))
  p2g_celltype_link[[i]] <- p2g_sub$Peak2GeneLinks
}

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'PTGER4', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(3, 10),
  upstream = 350000,
  downstream = 50000,
  loops = p2g_celltype_link,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label_link(p)
pdf('DOGMA_filtered_multiome/Plots/Peak2Gene-Track_modify_label_PTGER4_overlap_link_0.01_long.pdf', width = 10, height = 10)
p_new
dev.off()

############
p2g_celltype_link <- list()
for(i in names(table(proj$celltype_id))){
  p2g_sub <- readRDS(paste0('DOGMA_filtered_multiome/', i, '/p2g_fdr_0.001_loop.RDS'))
  p2g_celltype_link[[i]] <- p2g_sub$Peak2GeneLinks
}

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "celltype_id", 
  geneSymbol = 'PTGER4', 
  plotSummary = c("bulkTrack", "loopTrack"),
  sizes = c(3, 10),
  upstream = 350000,
  downstream = 50000,
  loops = p2g_celltype_link,
  useMatrix = 'GeneExpressionMatrix'
)

p_new <- modifyTrack_label_link(p)
pdf('DOGMA_filtered_multiome/Plots/Peak2Gene-Track_modify_label_PTGER4_overlap_link_0.001_long.pdf', width = 10, height = 10)
p_new
dev.off()
