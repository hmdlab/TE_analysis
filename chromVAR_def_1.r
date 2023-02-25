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
metadata <- read.table('data/metadata_cortex_blank_del.tsv',
                         header = TRUE,
                         stringsAsFactors=FALSE,quote="\t",row.names=1)#cell_label


peakfile<-"data/macs2/merged_cortex.bed"
peaks <- getPeaks(peakfile, sort_peaks = TRUE)

peaks <- resize(peaks, width = 500, fix = "center")

write.table(peaks,sep = "\t", "out/merged_cortex_500bp_pre.txt", quote=F, col.names=F, append=F)

bamfile <- "data/PreFrontalCortex_62216_RG_reheader.bam"

#gc(reset = TRUE)
#gc(reset = TRUE)
fragment_counts <- getCounts(bamfile, 
                             peaks, 
                             paired =  TRUE, 
                             by_rg = TRUE,
                             format = "bam",
                             colData = NULL)
name<-scan("data/metadata_cortex_label.txt", what = character(), sep = "\n", blank.lines.skip = FALSE, comment.char = "", quiet = TRUE)

fragment_counts_del<-fragment_counts[, fragment_counts@colData@rownames %in% name]
rm(fragment_counts)
gc(reset = TRUE)
gc(reset = TRUE)
filepath <- ".fa"
fasta <- readDNAStringSet(filepath, format = "fasta")

fragment_counts_del <- addGCBias(fragment_counts_del,genome=fasta)

counts_filtered <- filterPeaks(fragment_counts_del, non_overlapping = TRUE)
writeMM(counts_filtered@assays@data@listData[["counts"]],"out/counts_filtered.mtx")
write(counts_filtered@colData@rownames,"out/celllabel.txt")

bg <- getBackgroundPeaks(counts_filtered)

saveRDS(bg, file ='out/background_peaks_kmers_500bp.rds')#background


#kmer accessibility profile
kmer_ix_a <- matchKmers(7, counts_filtered,genome=fasta)
dev_a <- computeDeviations(object = counts_filtered, annotations = kmer_ix_a,
                         background_peaks = bg)

saveRDS(dev_a, file = 'out/Cortex_7mers_allgenome_zscore_500bp.rds')

#tSNE
dev_a@colData@rownames=sapply(strsplit(dev_a@colData@rownames,'.',fixed = TRUE), "[[", 1)
ix=match(dev_a@colData@rownames,rownames(metadata))
rownames(metadata)=rownames(metadata)[ix]
metadata$cell_label=metadata$cell_label[ix]

dev_a$label <- metadata$cell_label
tsne_results <- deviationsTsne(dev_a, threshold = 1.5, perplexity = 30)#parameter
tsne_plots <- plotDeviationsTsne(dev_a, tsne_results, 
                                 sample_column = "label", shiny = FALSE)
tsne_plots[[1]]
saveRDS(tsne_results, file = 'out/denovo_motifs_7mers_allgenome_tsne_results.rds')

kmername="CAGATGG"
tsne_plots <- plotDeviationsTsne(dev_a, tsne_results, annotation =kmername,
                                 sample_column = "label", shiny = FALSE)

tsne_plots[[1]]
tsne_plots[[2]]
ggsave(file = "out/tSNE_7mer_allgenome_30_denovo1.png", plot = tsne_plots[[2]], dpi = 400, width = 5.4, height = 3.6)

write.table(list(dev_a@colData@listData[["label"]])[[1]],"out/cell_type_label.txt")

#denovo motif
de_novos <- assembleKmers(dev_a, progress = FALSE,threshold = 2) #no progress bar
saveRDS(de_novos, file = 'out/denovo_motifs_7mers_allgenome.rds')



