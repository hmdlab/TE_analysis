#!/bin/bash
#$ -S /bin/bash
#$ -cwd

bedtools intersect -a out/denovo_1_motif_pos_ATAC.bed -b data/mm9.fa.align_map_onlyTE.bed -wa -wb > out/denovo_1_motif_pos_TEsubfamily_ATAC_align.bed
bedtools intersect -a out/denovo_3_motif_pos_ATAC.bed -b data/mm9.fa.align_map_onlyTE.bed -wa -wb > out/denovo_3_motif_pos_TEsubfamily_ATAC_align.bed

python3 Freq_score_TEdis_chromVAR.py -te MER130 -denovo 1
python3 Freq_score_TEdis_chromVAR.py -te MamRep434 -denovo 3
