bedtools intersect -a out/merged_cortex_500bp.txt -b data/data_3_bigdata_mm9_onlyTE.bed -wa -wb > out/merged_cortex_500bp_TE.bed
python3 TE_ATAC_rate.py
