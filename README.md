# TE_analysis
for "Transposons contribute to the acquisition of cell type-specific cis-elements in the brain"

[![DOI](https://zenodo.org/badge/606271716.svg)](https://zenodo.org/badge/latestdoi/606271716)

data source(scATAC-seq): https://atlas.gs.washington.edu/mouse-atac/

### Prepare
1. Get bam/sam file
```
$ cd data
$ wget http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/PreFrontalCortex_62216.bam
$ wget http://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/PreFrontalCortex_62216.bam.bai
$ samtools view -h data/PreFrontalCortex_62216.bam > data/PreFrontalCortex_62216.sam
```

2. Data shaping
```
# ATAC-seq
$ bash sam_RG_add.sh
$ bash sam_bunkatsu.sh
$ source ~/tools/MACS2/bin/activate   # MACS2 environment
$ macs2 callpeak -t data/PreFrontalCortex_62216_RG_reheader_Inhibitory_neurons.bam \
-f BAM --nomodel --keep-dup all --extsize 200 --shift -100 -g mm -n PreFrontalCortex_62216_RG_reheader_Inhibitory_neurons \
-B --outdir macs2
$ macs2 callpeak -t data/PreFrontalCortex_62216_RG_reheader_Ex._neurons_SCPN.bam \
-f BAM --nomodel --keep-dup all --extsize 200 --shift -100 -g mm -n PreFrontalCortex_62216_RG_reheader_Ex._neurons_SCPN \
-B --outdir macs2
$ macs2 callpeak -t data/PreFrontalCortex_62216_RG_reheader_Ex._neurons_CPN.bam \
-f BAM --nomodel --keep-dup all --extsize 200 --shift -100 -g mm -n PreFrontalCortex_62216_RG_reheader_Ex._neurons_CPN \
-B --outdir macs2
$ macs2 callpeak -t data/PreFrontalCortex_62216_RG_reheader_Ex._neurons_CThPN.bam \
-f BAM --nomodel --keep-dup all --extsize 200 --shift -100 -g mm -n PreFrontalCortex_62216_RG_reheader_Ex._neurons_CThPN \
-B --outdir macs2
$ macs2 callpeak -t data/PreFrontalCortex_62216_RG_reheader_SOM+_Interneurons.bam \
-f BAM --nomodel --keep-dup all --extsize 200 --shift -100 -g mm -n PreFrontalCortex_62216_RG_reheader_SOM+_Interneurons \
-B --outdir macs2
$ macs2 callpeak -t data/PreFrontalCortex_62216_RG_reheader_Microglia.bam \
-f BAM --nomodel --keep-dup all --extsize 200 --shift -100 -g mm -n PreFrontalCortex_62216_RG_reheader_Microglia \
-B --outdir macs2
$ macs2 callpeak -t data/PreFrontalCortex_62216_RG_reheader_Oligodendrocytes.bam \
-f BAM --nomodel --keep-dup all --extsize 200 --shift -100 -g mm -n PreFrontalCortex_62216_RG_reheader_Oligodendrocytes \
-B --outdir macs2
$ macs2 callpeak -t data/PreFrontalCortex_62216_RG_reheader_Astrocytes.bam \
-f BAM --nomodel --keep-dup all --extsize 200 --shift -100 -g mm -n PreFrontalCortex_62216_RG_reheader_Astrocytes \
-B --outdir macs2
$ deactivate
$ cat data/macs2/PreFrontalCortex_62216_RG_reheader_Astrocytes_peaks.narrowPeak data/macs2/PreFrontalCortex_62216_RG_reheader_Ex._neurons_CPN_peaks.narrowPeak data/macs2/PreFrontalCortex_62216_RG_reheader_Ex._neurons_CThPN_peaks.narrowPeak data/macs2/PreFrontalCortex_62216_RG_reheader_Ex._neurons_SCPN_peaks.narrowPeak data/macs2/PreFrontalCortex_62216_RG_reheader_Inhibitory_neurons_peaks.narrowPeak data/macs2/PreFrontalCortex_62216_RG_reheader_Microglia_peaks.narrowPeak data/macs2/PreFrontalCortex_62216_RG_reheader_Oligodendrocytes_peaks.narrowPeak data/macs2/PreFrontalCortex_62216_RG_reheader_SOM+_Interneurons_peaks.narrowPeak | sort -k1,1 -k2,2n | bedtools merge -i - > data/macs2/merged_cortex.bed

# TE(.out, .align)
$ cd data
$ wget https://www.repeatmasker.org/genomes/mm9/RepeatMasker-rm328-db20090604/mm9.fa.out.gz
$ wget https://www.repeatmasker.org/genomes/mm9/RepeatMasker-rm328-db20090604/mm9.fa.align.gz
$ cd ..
$ python3 TE_bedmake.py
```

※The format of data_3_bigdata_mm9_onlyTE.bed
| chr | start | end | TE subfamily | TE class|
| ------ | ---- | ---- | ---- | ---- |
| chr1	|3001723	|3002005	|RLTR25B	|LTR/ERVK|

※The format of mm9.fa.align_map_onlyTE.bed
| chr | start | end | TE subfamily | TE class| starting position of match in database sequence / (Complement) no. of bases in complement of the repeat consensus sequence prior to beginning of the match | ending position of match in database sequence / (Complement) starting position of match in database sequence | no. of bases in the repeat consensus sequence prior to beginning of the match / (Complement) ending position of match in database sequence |% of bases opposite a gap in the query sequence (deleted bp) | % of bases opposite a gap in the repeat consensus (inserted bp)|Line number of "mm9.fa.align"|
| ------ | ---- | ---- | ---- | ---- |---- | ---- | ---- |---- | ---- | ---- |
| chr1	|3001723	|3002005	|RLTR25B	|LTR/ERVK	|(0)	|1028	|625	|33.17	|0.72	|196|


※The format of (TF name)_peak_merge.bed
| chr | start | end | 
| ------ | ---- | ---- |
| chr1 | 100 | 200 | 

### De novo motifs with high variability in chromatin accessibility across cells are similar to known binding motifs of neural differentiation-related transcription factors.

- Visualization of cell similarity using t-SNE based on 7-mer chromatin accessibility with chromVAR (Figure 2a)
- Visualization of the accessibility of the seed k-mer used to generate the de novo motif with chromVAR (Figure 2b)
- Generating de novo motifs based on the k-mers with large accessibility variation across the cells with chromVAR (Figure 2c)
```
$ r chromVAR_def_1.r
```

-  Expression levels of scRNA-seq data based on the cell type labels transferred from scATAC-seq data with Seurat (Figure 2d)
```
$ cd  data
$ wget https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1
$ cd ..
$ r seurat_scRNA_scATAC.r
```

### Each Neurod2 and Lhx2-like motifs is accessible in a putative neuronal progenitor cell population in the adult brain.

Details can be found in [2.2.STREAM_scATAC-seq_k-mers.ipynb](https://github.com/pinellolab/STREAM/blob/master/tutorial/2.2.STREAM_scATAC-seq_k-mers.ipynb) (Figure 3) 

### Specific TE subfamilies are enriched in the accessible de novo motifs and transcription factor binding sites.

- Distribution of the percentage of ATAC peaks that overlap with TEs for each cell belonging to the cell type (Figure 4a,b)
```
$ bash TE_ATAC_rate.sh
```

- Distribution of classes of the TE that resides in transcription factor binding sites or ATAC-seq peaks (Figure 4c,d)
```
$ bash Distribution_of_classes_of_the_TE.sh
```

- Enrichment scores of the transcription factor ChIP peaks / Enrichment scores of the accessible motifs (Figure 5a) 
```
$ bash ESscore_TE_TF_neural_teleng.sh
```

- z-score using the enrichment score of accessible motifs
```
$ bash ESscore_TE_TF_neural_teleng_control.sh
```

-  TE deviation z-scores per cell based on the overlap of scATAC-seq peaks in the TE region with chromVAR (Figure 5b,c,d)
```
$ r chromVAR_TE.r
```

### Specific TEs, including MER130 and MamRep434, can function as important cis-elements of neuronal development genes.

- Distribution of transcription factor ChIP-seq reads mapped from the genome to consensus TE sequences (Figure 6b,c, S11) 
```
$ bash Freq_score_TEdis.sh
```

-  Distribution of the detection sites of accessible de novo motifs mapped from the genome to consensus TE sequences (Figure 6d, S11) 
```
$ bash Freq_score_TEdis_chromVAR.sh
$ Freq_score_TEdis_vis.py -tf (TF name) -te (TE name) -tfm (de novo motif number) # visualization
```


- IP/input ratio for the proportion of reads that mapped to TE-derived de novo motif sites (Figure S10)
```
$ bash deeptools.sh
```

### MER130 and MamRep434-derived cis-elements may contribute to brain evolution.

- Sequence conservation of the accessible motifs within TEs among mammals (Figure S13)
```
$ bash phyloP_socre_TE_motif.sh
```
