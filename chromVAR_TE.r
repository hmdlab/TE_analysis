#install.packages("BiocManager")
#BiocManager::install(c("Biostrings", "GenomicRanges", "rtracklayer"))
#if (!requireNamespace("BiocManager", quietly = TRUE))
 #   install.packages("BiocManager")

#BiocManager::install("TFBSTools")

#library(TFBSTools)

#if (!requireNamespace("BiocManager", quietly = TRUE))
 #   install.packages("BiocManager")

#BiocManager::install("chromVAR")


library("Biostrings")
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library('JASPAR2016')
library(BSgenome.Mmusculus.UCSC.mm9)
library(TFBSTools)
register(SerialParam())
#register(MulticoreParam(10))
Sys.setenv('R_MAX_VSIZE'=256000000000)
Sys.getenv('R_MAX_VSIZE')


counts_filtered<-readRDS("out/counts_filtered_peaks_kmers_500bp.rds")
bg<-readRDS("out/background_peaks_kmers_500bp.rds")

chip="data/data_3_bigdata_mm9_onlyTE.bed"
anno_ix <- getAnnotations(chip, 
                           rowRanges = rowRanges(counts_filtered), column = 4)
dev <- computeDeviations(object = counts_filtered, 
                                 annotations = anno_ix,background_peaks=bg)

saveRDS(dev, file = 'out/Cortex_TE_allgenome_zscore_500bp_onlycortex.rds')
