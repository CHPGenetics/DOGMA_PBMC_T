library(data.table)
library(parallel)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(venn)
library(ggvenn)
library(ggsci)
library(ggrepel)
library(metafor)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/code')

load('../plots/eQTA/pseudobulk_eQTA_re.RData')
eqta <- me$cis$eqtls
eqta$pair <- paste(eqta$gene, eqta$snps, sep = ':')

load('~/RWorkSpace/Duerr_bulk/eQTA/eQTA.RData')
bulk <- me$cis$eqtls
bulk$pair <- paste(bulk$gene, bulk$snps, sep = ':')

draw_plot <- function(df1, df2){
  df <- merge(df1, df2, by = 'pair')
  df <- df[order(df$pvalue.x),]
  df$sign <- 'Not Significant'
  df$sign[df$FDR.x < 0.05] <- 'Significant'
  df$sign[df$FDR.x < 0.05 & df$FDR.y < 0.05 & df$beta.x*df$beta.y > 0] <- 'Replicated'
  df$sign <- factor(df$sign, levels = c('Replicated', 'Significant', 'Not Significant'))
  p <- ggplot(df, aes(y = beta.x, x = beta.y))+
    geom_point(aes(color = sign))+
    scale_color_manual(values = c("firebrick", alpha("dodgerblue3", 0.7), alpha("lightgrey", 0.1))) +
    theme_bw(base_size = 12) +
    theme(text = element_text(size=16)) +
    theme(legend.position = "bottom") + 
    xlab('Single-cell coefficient')+
    ylab('Bulk coefficient')+
    labs(col = '', subtitle = paste0(table(df$sign)['Replicated'], ' out of ', table(df$sign)['Replicated']+table(df$sign)['Significant'], ' significant eQTA pairs were replicated in bulk'))+
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

p <- draw_plot(eqta, bulk)

pdf('../plots/pseudo_bulk/merge/eQTA/4.pdf', width = 7, height = 7)
p
dev.off()
