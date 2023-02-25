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
register(SerialParam())
#register(MulticoreParam(10))
Sys.setenv('R_MAX_VSIZE'=256000000000)
Sys.getenv('R_MAX_VSIZE')

library(Biostrings)
library(ggmotif) 
library(TFBSTools)

PWN<-readRDS("out/denovo_motifs_7mers_allgenome.rds")
## motif scan
filepath <- ".fa"
fasta <- readDNAStringSet(filepath, format = "fasta")
motif_pos <- matchMotifs(PWM,counts_filtered_peaks_kmers_500bp@rowRanges,genome=fasta,
                         out = "positions") 

write.table(motif_pos$denovo_1, "out/denovo_1_motif_pos.txt", quote=F, 
                                  col.names=T, append=F)
write.table(motif_pos$denovo_2, "out/denovo_2_motif_pos.txt", quote=F, 
                                  col.names=T, append=F)
write.table(motif_pos$denovo_3, "out/denovo_3_motif_pos.txt", quote=F, 
                                  col.names=T, append=F)
write.table(motif_pos$denovo_4, "out/denovo_4_motif_pos.txt", quote=F, 
                                  col.names=T, append=F)
