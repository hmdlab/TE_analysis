#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#calculate
python3 ESscore_TE_TF_neural_teleng_control.py -tf denovo_1 -pm m -pre Y
for i in `seq 50`
do
    bedtools intersect -a out/denovo_1_motif_pos_ATAC_merge_control_$i.bed -b data/data_3_bigdata_mm9_onlyTE.bed -wa -wb > out/denovo_1_motif_pos_TEsubfamily_ATA_control_$i.txt
done
python3 ESscore_TE_TF_neural_teleng_control.py -tf denovo_1 -pm m

python3 ESscore_TE_TF_neural_teleng_control.py -tf denovo_2 -pm m -pre Y
for i in `seq 50`
do
    bedtools intersect -a out/denovo_2_motif_pos_ATAC_merge_control_$i.bed -b data/data_3_bigdata_mm9_onlyTE.bed -wa -wb >out/denovo_2_motif_pos_TEsubfamily_ATA_control_$i.txt
done
python3 ESscore_TE_TF_neural_teleng_control.py -tf denovo_2 -pm m

python3 ESscore_TE_TF_neural_teleng_control.py -tf denovo_3 -pm m -pre Y
for i in `seq 50`
do
    bedtools intersect -a out/denovo_3_motif_pos_ATAC_merge_control_$i.bed -b data/data_3_bigdata_mm9_onlyTE.bed -wa -wb > out/denovo_3_motif_pos_TEsubfamily_ATA_control_$i.txt
done
python3 ESscore_TE_TF_neural_teleng_control.py -tf denovo_3 -pm m

python3 zscore.py -denovo 1 -pse 0.5
python3 zscore.py -denovo 2 -pse 0.5
python3 zscore.py -denovo 3 -pse 0.5
