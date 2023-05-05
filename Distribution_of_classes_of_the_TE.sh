#!/bin/bash
#$ -S /bin/bash
#$ -cwd

# prepare
r chromVAR_motif_scan.r
python3 downstream.py
cat out/merged_cortex_500bp_pre.txt  | cut -f 2,3,4 > out/merged_cortex_500bp.txt
bedtools intersect -a out/denovo_1_motif_pos.txt -b out/merged_cortex_500bp.txt -wa -wb > out/denovo_1_motif_pos_ATAC_merge.bed
bedtools intersect -a out/denovo_2_motif_pos.txt -b out/merged_cortex_500bp.txt -wa -wb > out/denovo_2_motif_pos_ATAC_merge.bed
bedtools intersect -a out/denovo_3_motif_pos.txt -b out/merged_cortex_500bp.txt -wa -wb > out/denovo_3_motif_pos_ATAC_merge.bed
bedtools intersect -a out/denovo_1_motif_pos_ATAC_merge.bed -b data/data_3_bigdata_mm9_onlyTE.bed -wa -wb > out/denovo_1_motif_pos_TEsubfamily_ATAC.bed
bedtools intersect -a out/denovo_2_motif_pos_ATAC_merge.bed -b data/data_3_bigdata_mm9_onlyTE.bed -wa -wb > out/denovo_2_motif_pos_TEsubfamily_ATAC.bed
bedtools intersect -a out/denovo_3_motif_pos_ATAC_merge.bed -b data/data_3_bigdata_mm9_onlyTE.bed -wa -wb > out/denovo_3_motif_pos_TEsubfamily_ATAC.bed


python3 Distribution_of_classes_of_the_TE.py
