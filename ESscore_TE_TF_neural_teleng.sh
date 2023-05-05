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


# calculate
bedtools intersect -a data/data_3_bigdata_mm9_onlyTE.bed -b data/TF_chipatlas_Neurod2_merge.bed -wa -wb > data/TE_contain_Neurod2_peak_neural_merge.bed
python3 ESscore_TE_TF_neural_teleng.py -tf Neurod2 -pm p
python3 ESscore_TE_TF_neural_teleng.py -tf denovo_1 -pm m

python3 ESscore_TE_TF_neural_teleng.py -tf denovo_2 -pm m

bedtools intersect -a data/data_3_bigdata_mm9_onlyTE.bed -b data/TF_chipatlas_Lhx2_merge.bed -wa -wb > data/TE_contain_Lhx2_peak_neural_merge.bed
python3 ESscore_TE_TF_neural_teleng.py -tf Lhx2 -pm p
python3 ESscore_TE_TF_neural_teleng.py -tf denovo_3 -pm m
