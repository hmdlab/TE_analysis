#!/bin/bash
#$ -S /bin/bash
#$ -cwd


python3 sam_RG_add.py
samtools view -Sb data/PreFrontalCortex_62216_RG.sam > data/PreFrontalCortex_62216_RG.bam
samtools index data/PreFrontalCortex_62216_RG.bam
samtools reheader data/PreFrontalCortex_62216_RG_header.sam data/PreFrontalCortex_62216_RG.bam > data/PreFrontalCortex_62216_RG_reheader.bam
samtools index data/PreFrontalCortex_62216_RG_reheader.bam

sinto fragments -b data/PreFrontalCortex_62216_RG_reheader.bam -f data/PreFrontalCortex_62216_RG_reheader_fragments.bed -t RG
sort -k1,1 -k2,2n data/PreFrontalCortex_62216_RG_reheader_fragments.bed > data/PreFrontalCortex_62216_RG_reheader_fragments.sort.bed
bgzip data/PreFrontalCortex_62216_RG_reheader_fragments.sort.bed
tabix -p bed data/PreFrontalCortex_62216_RG_reheader_fragments.sort.bed.gz

