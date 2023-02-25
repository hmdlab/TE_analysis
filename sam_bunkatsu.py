#!/usr/bin/env python
# coding: utf-8

import csv
import numpy as np


import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-celltype','--celltype',type=str,default='NA',help='celltype')
args=parser.parse_args()
celltype=args.celltype

import pysam
ATAC = pysam.AlignmentFile("data/PreFrontalCortex_62216_RG_reheader.sam", "r")

label=[]
#for read in ATAC:
 #   print(str(read).split("\t"))
  #  break
for read in ATAC:
    label.append(str(read).split("\t")[0].split(':')[0])

#print(len(label))
ct=len(label)
label_name=list(set(label))
#print(len(label),len(label_name),label_name[0],label[0])


with open('data/metadata_cortex_blank_del.tsv') as f:
    reader = csv.reader(f, delimiter='\t')
    cell = [row for row in reader]
cellid=[]
celltypeid=[]
for i in range(1,len(cell)):
  cellid.append(cell[i][0])
  celltypeid.append(cell[i][1])


in_samfile=pysam.AlignmentFile("data/PreFrontalCortex_62216_RG_reheader.sam", "r")
out_samfile = pysam.AlignmentFile("data/PreFrontalCortex_62216_RG_reheader_"+str(celltype)+".sam", "w", template = in_samfile)
j=0
for j in range(len(label)):
    read=next(in_samfile)
    if cellid.count(str(read).split("\t")[0].split(':')[0])>0:
        if celltypeid[cellid.index(str(read).split("\t")[0].split(':')[0])]==celltype:
            out_samfile.write(read)
    j+=1
    if j==ct:
        break
in_samfile.close()
out_samfile.close()
