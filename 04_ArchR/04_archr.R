library(ArchR)
library(stringr)
set.seed(1)

setwd('~/RWorkSpace/DOGMA-seq/PBMC/output/ArchR/')

addArchRThreads(threads = 64) 

paths <- c('Duerr_20210419_DOGMAseq_DIG_arc',
           'Duerr_20210610_DOGMAseq-1_arc', 'Duerr_20210610_DOGMAseq-2_arc',
           'Duerr_20210826_DOGMAseq-1_arc', 'Duerr_20210826_DOGMAseq-2_arc',
           'Duerr_20210831_DOGMAseq-1_arc', 'Duerr_20210831_DOGMAseq-2_arc',
           'Duerr_20220126_DOGMAseq_arc',
           'Duerr_20220215_DOGMAseq-1_arc', 'Duerr_20220215_DOGMAseq-2_arc',
           'Duerr_20220218_DOGMAseq-1_arc', 'Duerr_20220218_DOGMAseq-2_arc',
           'Duerr_20220524_DOGMAseq-1_arc', 'Duerr_20220524_DOGMAseq-2_arc',
           'Duerr_20220603_DOGMAseq-2_arc'
)

inputFiles <- paste0('../../../', paths, '/atac_fragments.tsv.gz')

names(inputFiles) <- c('04191',
                       '06101', '06102',
                       '08261', '08262',
                       '08311', '08312',
                       '01261',
                       '02151', '02152',
                       '02181', '02182',
                       '05241', '05242',
                       '06032')

addArchRGenome("hg38")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 0, #Dont set this too high because you can always increase later
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = 'DOGMA',
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA')

meta <- read.csv('../tcell_annotated_updated.csv', row.names = 'X')

meta$cellNames <- rownames(meta)
meta$cellNames <- gsub('_', '#', meta$cellNames)

table(proj$cellNames %in% meta$cellNames)
proj <- proj[proj$cellNames %in% meta$cellNames,]
rownames(meta) <- meta$cellNames
meta <- meta[rownames(proj),]
proj$condition <- meta$condition
proj$sample <- meta$sample
proj$celltype <- meta$celltype_updated

doubScores <- addDoubletScores(
  input = proj,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = 'DOGMA_filtered')
