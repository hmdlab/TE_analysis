#!/bin/bash
#$ -S /bin/bash
#$ -cwd


# 1. get fastq
#Neurod2
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR3955800/SRR3955800
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5170934/SRR5170934
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5170935/SRR5170935
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5903318/SRR5903318

#Lhx2
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR3955796/SRR3955796
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR3955797/SRR3955797
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR3955798/SRR3955798
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5170930/SRR5170930
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5170931/SRR5170931
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5170940/SRR5170940
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5903317/SRR5903317


fastq-dump SRR3955800
fastq-dump SRR5170934 --split-files
fastq-dump SRR5170935 --split-files
fastq-dump SRR5903318


fastq-dump SRR3955796
fastq-dump SRR3955797
fastq-dump SRR3955798
fastq-dump SRR5170930 --split-files
fastq-dump SRR5170931 --split-files
fastq-dump SRR5170940
fastq-dump SRR5903317


#rm SRR3955796
#rm SRR3955797
#rm SRR3955798

#rm SRR5170930
#rm SRR5170931
#rm SRR5170940
#rm SRR5903317
#rm SRR3955800
#rm SRR5170934
#rm SRR5170935
#rm SRR5903318

# 2. mappping
mkdir bowtie2_index    # Pre-built index
cd bowtie2_index   
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm9.zip
unzip mm9.zip
cd ..
bowtie2 -p 2 -x bowtie2_index/mm9 \
    -U SRR3955796.fastq | samtools view -F 0x04 -Sb - > SRR3955796.bam

samtools sort SRR3955796.bam -o SRR3955796.sorted.bam
#rm SRR3955796.fastq
#rm SRR3955796.bam

bowtie2 -p 2 -x bowtie2_index/mm9 \
    -U SRR3955797.fastq | samtools view -F 0x04 -Sb - > SRR3955797.bam

samtools sort SRR3955797.bam -o SRR3955797.sorted.bam
#rm SRR3955797.fastq
#rm SRR3955797.bam

bowtie2 -p 2 -x bowtie2_index/mm9 \
    -U SRR3955798.fastq | samtools view -F 0x04 -Sb - > SRR3955798.bam

samtools sort SRR3955798.bam -o SRR3955798.sorted.bam
#rm SRR3955798.fastq
#rm SRR3955798.bam


bowtie2 -p 2 -x bowtie2_index/mm9 \
    -1 SRR5170930_1.fastq -2 SRR5170930_2.fastq | samtools view  -F 0x04 -Sb - > SRR5170930.bam

samtools sort SRR5170930.bam -o SRR5170930.sorted.bam
#rm SRR5170930_1.fastq
#rm SRR5170930_2.fastq
#rm SRR5170930.bam

bowtie2 -p 2 -x bowtie2_index/mm9 \
    -1 SRR5170931_1.fastq -2 SRR5170931_2.fastq | samtools view  -F 0x04 -Sb - > SRR5170931.bam

samtools sort SRR5170931.bam -o SRR5170931.sorted.bam
#rm SRR5170931_1.fastq
#rm SRR5170931_2.fastq
#rm SRR5170931.bam

bowtie2 -p 2 -x bowtie2_index/mm9 \
    -U SRR5170940.fastq | samtools view -F 0x04 -Sb - > SRR5170940.bam

samtools sort SRR5170940.bam -o SRR5170940.sorted.bam
#rm SRR5170940.fastq
#rm SRR5170940.bam

bowtie2 -p 2 -x bowtie2_index/mm9 \
    -U SRR5903317.fastq | samtools view -F 0x04 -Sb - > SRR5903317.bam

samtools sort SRR5903317.bam -o SRR5903317.sorted.bam
#rm SRR5903317.fastq
#rm SRR5903317.bam

bowtie2 -p 2 -x bowtie2_index/mm9 \
    -1 SRR5170935_1.fastq -2 SRR5170935_2.fastq | samtools view  -F 0x04 -Sb - > SRR5170935.bam

samtools sort SRR5170935.bam -o SRR5170935.sorted.bam
#rm SRR5170935_1.fastq
#rm SRR5170935_2.fastq
#rm SRR5170935.bam

bowtie2 -p 2 -x bowtie2_index/mm9 \
    -U SRR5903318.fastq | samtools view  -F 0x04 -Sb - > SRR5903318.bam

samtools sort SRR5903318.bam -o SRR5903318.sorted.bam
#rm SRR5903318.fastq
#rm SRR5903318.bam


samtools rmdup -s SRR3955800.sorted.bam SRR3955800.sorted.rmdup.bam
samtools rmdup SRR5170934.sorted.bam SRR5170934.sorted.rmdup.bam
samtools rmdup SRR5170935.sorted.bam SRR5170935.sorted.rmdup.bam
samtools rmdup -s SRR5903318.sorted.bam SRR5903318.sorted.rmdup.bam  
samtools rmdup -s SRR3955796.sorted.bam SRR3955796.sorted.rmdup.bam
samtools rmdup -s SRR3955797.sorted.bam SRR3955797.sorted.rmdup.bam
samtools rmdup -s SRR3955798.sorted.bam SRR3955798.sorted.rmdup.bam

samtools rmdup -s SRR5170930.sorted.bam SRR5170930.sorted.rmdup.bam
samtools rmdup -s SRR5170931.sorted.bam SRR5170931.sorted.rmdup.bam
samtools rmdup -s SRR5170940.sorted.bam SRR5170940.sorted.rmdup.bam
samtools rmdup -s SRR5903317.sorted.bam SRR5903317.sorted.rmdup.bam

samtools index data/SRR3955800.sorted.rmdup.bam 
samtools index data/SRR3955796.sorted.rmdup.bam 
samtools index data/SRR3955797.sorted.rmdup.bam 
samtools index data/SRR3955798.sorted.rmdup.bam 
samtools index data/SRR5170934.sorted.rmdup.bam 
samtools index data/SRR5170935.sorted.rmdup.bam 
samtools index data/SRR5903318.sorted.rmdup.bam 
samtools index data/SRR5170930.sorted.rmdup.bam 
samtools index data/SRR5170931.sorted.rmdup.bam 
samtools index data/SRR5170940.sorted.rmdup.bam 
samtools index data/SRR5903317.sorted.rmdup.bam


bedtools intersect -a out/denovo_1_motif_pos.txt -b data/data_3_bigdata_mm9_onlyTE.bed -wa -wb > out/denovo_1_motif_pos_TEsubfamily.txt
bedtools intersect -a out/denovo_3_motif_pos.txt -b data/data_3_bigdata_mm9_onlyTE.bed -wa -wb > out/denovo_3_motif_pos_TEsubfamily.txt

python3 deeptools_bedmake.py


module load singularity
cd /usr/local/biotools/d
singularity exec deeptools:3.5.1--py_0 plotEnrichment -b data/SRR3955800.sorted.rmdup.bam data/SRR3955796.sorted.rmdup.bam data/SRR3955797.sorted.rmdup.bam data/SRR3955798.sorted.rmdup.bam \
--BED out/denovo_1_motif_pos_MER130.bed \
--regionLabels "de novo motif 1 in MER130" --variableScales --labels "input" "ip_1" "ip_2" "ip_3" --colors black black black black \
-o out/MER130_denovo1_enrichment.png  --outRawCounts out/MER130_denovo1_enrichment.txt

singularity exec deeptools:3.5.1--py_0 plotEnrichment -b data/SRR5170934.sorted.rmdup.bam data/SRR5170935.sorted.rmdup.bam data/SRR5903318.sorted.rmdup.bam data/SRR5170930.sorted.rmdup.bam data/SRR5170931.sorted.rmdup.bam data/SRR5170940.sorted.rmdup.bam data/SRR5903317.sorted.rmdup.bam \
--BED out/denovo_3_motif_pos_MamRep434.bed \
--regionLabels "de novo motif 3 in MamRep434" --variableScales --labels "input_1" "input_2" "input_3" "ip_1" "ip_2" "ip_3(KO)" "ip_4(KO)" --colors black black black black black black black \
-o out/MamRep434_denovo3_enrichment.png  --outRawCounts out/MamRep434_denovo3_enrichment.txt