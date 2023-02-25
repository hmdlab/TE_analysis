#!/bin/bash
#$ -S /bin/bash
#$ -cwd

cd data
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y
do
  wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/phyloP30way/placentalMammals/chr$i.phyloP30way.placental.wigFix.gz
done
cd ..

bedtools intersect -a data/MER130_region.bed -b data/macs2/merged_cortex.bed > out/MER130_ATAC_intersect.bed
bedtools intersect -a data/MamRep434_region.bed -b data/macs2/merged_cortex.bed > out/MamRep434_ATAC_intersect.bed
bedtools shuffle -i data/control_1000_7bp.bed -g data/mm9.genome.txt -noOverlapping -incl out/MER130_ATAC_intersect.bed > out/MER130_ATAC_intersect_control.bed
bedtools shuffle -i data/control_1000_7bp.bed -g data/mm9.genome.txt -noOverlapping -incl out/MamRep434_ATAC_intersect.bed > out/MamRep434_ATAC_intersect_control.bed
python3 phyloP_score_TE_motif.py -tf denovo_1 -te MER130
python3 phyloP_score_TE_motif.py -tf denovo_3 -te MamRep434
python3 phyloP_score_TE_motif.py -tf denovo_1 -te MER130 -cont True
python3 phyloP_score_TE_motif.py -tf denovo_3 -te MamRep434 -cont True
