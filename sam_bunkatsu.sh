#!/bin/bash
#$ -S /bin/bash
#$ -cwd

python3 sam_bunkatsu.py -celltype Inhibitory_neurons
python3 sam_bunkatsu.py -celltype Ex._neurons_SCPN
python3 sam_bunkatsu.py -celltype Ex._neurons_CPN
python3 sam_bunkatsu.py -celltype Ex._neurons_CThPN
python3 sam_bunkatsu.py -celltype SOM+_Interneurons
python3 sam_bunkatsu.py -celltype Microglia
python3 sam_bunkatsu.py -celltype Oligodendrocytes
python3 sam_bunkatsu.py -celltype Astrocytes
samtools view -Sb data/PreFrontalCortex_62216_RG_reheader_Inhibitory_neurons.sam > data/PreFrontalCortex_62216_RG_reheader_Inhibitory_neurons.bam
samtools view -Sb data/PreFrontalCortex_62216_RG_reheader_Ex._neurons_SCPN.sam > data/PreFrontalCortex_62216_RG_reheader_Ex._neurons_SCPN.bam
samtools view -Sb data/PreFrontalCortex_62216_RG_reheader_Ex._neurons_CPN.sam > data/PreFrontalCortex_62216_RG_reheader_Ex._neurons_CPN.bam
samtools view -Sb data/PreFrontalCortex_62216_RG_reheader_Ex._neurons_CThPN.sam > data/PreFrontalCortex_62216_RG_reheader_Ex._neurons_CThPN.bam
samtools view -Sb data/PreFrontalCortex_62216_RG_reheader_SOM+_Interneurons.sam > data/PreFrontalCortex_62216_RG_reheader_SOM+_Interneurons.bam
samtools view -Sb data/PreFrontalCortex_62216_RG_reheader_Microglia.sam > data/PreFrontalCortex_62216_RG_reheader_Microglia.bam
samtools view -Sb data/PreFrontalCortex_62216_RG_reheader_Oligodendrocytes.sam > data/PreFrontalCortex_62216_RG_reheader_Oligodendrocytes.bam
samtools view -Sb data/PreFrontalCortex_62216_RG_reheader_Astrocytes.sam > data/PreFrontalCortex_62216_RG_reheader_Astrocytes.bam
samtools index data/PreFrontalCortex_62216_RG_reheader_Inhibitory_neurons.bam
samtools index data/PreFrontalCortex_62216_RG_reheader_Ex._neurons_SCPN.bam
samtools index data/PreFrontalCortex_62216_RG_reheader_Ex._neurons_CPN.bam
samtools index data/PreFrontalCortex_62216_RG_reheader_Ex._neurons_CThPN.bam
samtools index data/PreFrontalCortex_62216_RG_reheader_SOM+_Interneurons.bam
samtools index data/PreFrontalCortex_62216_RG_reheader_Microglia.bam
samtools index data/PreFrontalCortex_62216_RG_reheader_Oligodendrocytes.bam
samtools index data/PreFrontalCortex_62216_RG_reheader_Astrocytes.bam
