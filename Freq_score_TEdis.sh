#!/bin/bash
#$ -S /bin/bash
#$ -cwd


bedtools intersect -a data/mm9.fa.align_map_onlyTE.bed -b data/Neurod2_peak_merge.bed -wa -wb > out/Neurod2_peak_mm9.fa.align_map_only_TE.bed
bedtools intersect -a data/mm9.fa.align_map_onlyTE.bed -b data/Lhx2_peak_merge.bed -wa -wb > out/Lhx2_peak_mm9.fa.align_map_only_TE.bed

python3 Freq_score_TEdis.py -te MER130 -tf Neurod2 -mask N
python3 Freq_score_TEdis.py -te MamRep434 -tf Lhx2 -mask N

