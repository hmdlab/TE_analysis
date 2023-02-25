#!/usr/bin/env python
# coding: utf-8

import csv
import numpy as np


import pysam


with open('data/metadata_cortex_blank_del.tsv') as f:
    reader = csv.reader(f, delimiter='\t')
    cell = [row for row in reader]

cellid=[]
celltypeid=[]
for i in range(1,len(cell)):
  cellid.append(cell[i][0])
  celltypeid.append(cell[i][1])

test = pysam.AlignmentFile("data/PreFrontalCortex_62216.sam", "r")
out_samfile = pysam.AlignmentFile("data/PreFrontalCortex_62216_RG.sam", "w", template = test)
for count, read in enumerate(test.fetch()):
    if cellid.count(str(read).split("\t")[0].split(':')[0])>0:
        new_tags = read.tags
        new_tags.append(('RG', cellid[cellid.index(str(read).split("\t")[0].split(':')[0])]))
        read.tags = new_tags
        out_samfile.write(read)
out_samfile.close()
