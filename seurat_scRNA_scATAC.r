library(SeuratData)

library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)

#for github
metadata <- read.table('data/metadata_cortex_blank_del.tsv',
                       header = TRUE,
                       stringsAsFactors=FALSE,quote="\t",row.names=1)#cell_label

#k-mer(scATAC-seq)
counts_filtered_peaks_kmers_500bp<-readRDS("out/Cortex_7mers_allgenome_zscore_500bp.rds")
# RDS(Allen,scRNA-seq)
pbmc.rna <- readRDS("data/allen_cortex.rds")

counts_filtered_peaks_kmers_500bp@colData@rownames<-sapply(strsplit(counts_filtered_peaks_kmers_500bp@colData@rownames,'.',fixed = TRUE), "[[", 1)
counts_filtered_peaks_kmers_500bp@assays@data@listData[["counts"]]@Dimnames[[2]]<-sapply(strsplit(counts_filtered_peaks_kmers_500bp@assays@data@listData[["counts"]]@Dimnames[[2]],'.',fixed = TRUE), "[[", 1)
ix=match(counts_filtered_peaks_kmers_500bp@colData@rownames,rownames(metadata))
rownames(metadata)=rownames(metadata)[ix]
metadata$cell_label=metadata$cell_label[ix]

chromatinassay <- CreateChromatinAssay(counts = counts_filtered_peaks_kmers_500bp@assays@data@listData[["counts"]], 
                                       genome = "mm9", ranges =counts_filtered_peaks_kmers_500bp@rowRanges)
pbmc.atac <- CreateSeuratObject(
  counts = chromatinassay,
  assay = "peaks",
)

# Perform standard analysis of each modality independently RNA analysis
pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna)
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna)
pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)

gtf <- rtracklayer::import('data/Mus_musculus.NCBIM37.67.gtf')
seqlevelsStyle(gtf) <- 'UCSC'
Annotation(pbmc.atac) <- gtf

fragments <- CreateFragmentObject(
    path = "data/PreFrontalCortex_62216_RG_reheader_fragments.sort.bed.gz",
    cells = colnames(pbmc.atac),
    validate.fragments = FALSE
  )

Fragments(pbmc.atac) <- fragments
activity.matrix <- GeneActivity(pbmc.atac, 
                               extend.upstream = 2000, verbose = TRUE)#No fragment information found for requested assay
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)

# We exclude the first dimension as this is typically correlated with sequencing depth
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")
pbmc.atac <- RunSVD(pbmc.atac)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
pbmc.atac$cell.type.label<-metadata$cell_label

p1 <- DimPlot(pbmc.rna, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(pbmc.atac, group.by = "cell.type.label", label = TRUE) + NoLegend() + ggtitle("ATAC")
p1 + p2


# normalize gene activities
DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac, features = rownames(pbmc.atac))

# Identify anchors
#ATAC→RNA
pbmc.atac <- FindVariableFeatures(pbmc.atac)
transfer.anchors <- FindTransferAnchors(reference = pbmc.atac, query = pbmc.rna, features = VariableFeatures(object = pbmc.atac),
                                        reference.assay = "ACTIVITY",query.assay = "RNA", reduction = "cca")#reference.assay = "RNA"


predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc.atac$cell.type.label,
  weight.reduction = pbmc.rna[['pca']],
  dims = 1:30
)

pbmc.rna <- AddMetaData(object = pbmc.rna, metadata = predicted.labels)


plot1 <- DimPlot(
  object = pbmc.atac,
  group.by = 'cell.type.label',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot2 <- DimPlot(
  object = pbmc.rna,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot1 + plot2



writeMM(pbmc.rna@assays[["RNA"]]@counts,"out/gene_expression_transfer_counts.mtx")
writeMM(pbmc.rna@assays[["RNA"]]@data,"out/gene_expression_transfer.mtx")
write.table(pbmc.rna@meta.data[["predicted.id"]],"out/gene_expression_transfer_celltypeid.txt")
write.table(pbmc.rna@assays[["RNA"]]@counts@Dimnames[[1]],"out/gene_expression_transfer_genenaame.txt")
write.table(pbmc.rna@assays[["gene"]]@counts@Dimnames[[2]],"/out/gene_expression_transfer_cellname.txt")

