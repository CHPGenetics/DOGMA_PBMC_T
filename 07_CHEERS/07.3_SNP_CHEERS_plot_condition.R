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

read.dir <- function(path, pattern){
  files <- list.files(path, pattern)
  tables <- list()
  for (i in files){
    tmp <- read.table(paste(path, i, sep = '/'), sep = '\t')
    colnames(tmp) <- c('name', 'p_value')
    tmp$SNP <- gsub(pattern, '', i)
    tables[[gsub(pattern, '', i)]] <- tmp
  }
  return(tables)
}

table <- read.dir('../output/SNP_2/CHEERS/AllPeaks/Merged_condition_analysis/', '_disease_enrichment_pValues.txt')

data_plot_bp <- table[[1]]
for(i in 2:length(table)){
  data_plot_bp <- rbind(data_plot_bp, table[[i]])
}

############
data_plot <- data_plot_bp[grep('IBD|UC|CD', data_plot_bp$SNP),]
data_plot$snp <- 'All variants'
data_plot$snp[grep('Lead', data_plot$SNP)] <- 'Leading variants'
data_plot$snp[grep('Causal', data_plot$SNP)] <- 'Causal variants with > 95% certainty'
data_plot$snp <- factor(data_plot$snp, levels = c('Causal variants with > 95% certainty', 'Leading variants', 'All variants'))
data_plot$group <- 'IBD'
data_plot$group[grep('UC', data_plot$SNP)] <- 'UC'
data_plot$group[grep('CD', data_plot$SNP)] <- 'CD'
data_plot$snp_group <- paste0(data_plot$snp, ' (', data_plot$group, ')')

data_plot$condition <- str_split_fixed(data_plot$name, '_', 3)[,3]
data_plot$cell <- paste(str_split_fixed(data_plot$name, '_', 3)[,1], str_split_fixed(data_plot$name, '_', 3)[,2], sep = '_')

###########
data_plot$snp <- gsub('Leading', 'Lead', data_plot$snp)
data_plot$snp_group <- gsub('Leading', 'Lead', data_plot$snp_group)
data_plot <- data_plot[data_plot$snp != 'Causal variants with > 95% certainty',]
data_plot$condition <- factor(data_plot$condition, levels = c('Act_IL1B_IL23', 'Act_IL1B_IL23_PGE2'))
data_plot$snp_group <- factor(data_plot$snp_group, levels = c('All variants (CD)', 'All variants (UC)', 'All variants (IBD)', 
                                                              'Lead variants (CD)', 'Lead variants (UC)', 'Lead variants (IBD)'))

data_plot$log10_p_value <- -log10(data_plot$p_value)
data_plot$log10_p_value[data_plot$p_value >= 0.05] <- NA

p <- ggplot(data_plot, 
            aes(y = condition, x = cell, fill = log10_p_value))+
  geom_tile()+
  scale_fill_viridis(option = "C", direction = -1, na.value = 'darkgrey', limits = c(1,4), oob = scales::squish)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'right', plot.title = element_text(hjust = 0.5))+
  labs(x = '', y = '', fill = '-log10(P-value)')+ 
  coord_fixed(expand = T)
p <- facet(p, facet.by = 'snp_group', ncol = 1, strip.position = 'left')+
  theme(strip.text.y.left = element_text(angle = 0), strip.background = element_blank())

pdf('../plots/SNP_2/CHEERS/enrich_disease_IBD_updated.pdf', width = 10, height = 7)
p
dev.off()

p <- facet(p, facet.by = 'snp_group', ncol = 1, strip.position = 'left')+
  theme(strip.text.y.left = element_text(angle = 0), strip.background = element_blank())

pdf('../plots/SNP_2/CHEERS/enrich_disease_IBD_updated_no_legend.pdf', width = 8, height = 6.5)
p+theme(legend.position = 'none')
dev.off()
