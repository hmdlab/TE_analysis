#!/usr/bin/env python
# coding: utf-8

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-denovo','--denovo',type=str,default='None',help='denovo motif_number')
parser.add_argument('-pse','--pse',type=float,default=0.5,help='pseudocount')
args=parser.parse_args()
motif=args.denovo
pse=args.pse

import math


f = open('data/te_tf_mat_TElabel_cortex.txt')
data3 = f.read() 
f.close()
tename= data3.split('\n')





import csv
import numpy as np


with open('out/denovo_'+str(motif)+'_TE_log2fold_count_motif_detail.bed') as f:
    reader = csv.reader(f, delimiter='\t')
    enrich = [row for row in reader]
enrich=np.array(enrich)
enrich=enrich.astype('float32')
enrich=enrich.T




score=[]
for i in range(len(tename)):
    if pse>0:
        score.append(math.log2((enrich[i][0]/enrich[i][1])/(enrich[i][2]/enrich[i][3])+pse))
    else:
        if enrich[i][0]!=0:#no pseudocount
            score.append(math.log2((enrich[i][0]/enrich[i][1])/(enrich[i][2]/enrich[i][3])))
        else:
            score.append(float('nan'))

score_TE=score


control=[]
for j in range(1,51):
    with open('out/denovo_'+str(motif)+'_TE_log2fold_count_motif_control_'+str(j)+'_detail.bed') as f:
        reader = csv.reader(f, delimiter='\t')
        enrich = [row for row in reader]
    enrich=np.array(enrich)
    enrich=enrich.astype('float32')
    enrich=enrich.T
    score=[]
    for i in range(len(tename)):
        if pse>0:
            score.append(math.log2((enrich[i][0]/enrich[i][1])/(enrich[i][2]/enrich[i][3])+pse))
        else:
            if enrich[i][0]!=0:#no pseudocount
                score.append(math.log2((enrich[i][0]/enrich[i][1])/(enrich[i][2]/enrich[i][3])))
            else:
                score.append(float('nan'))
    control.append(score)

control=np.array(control)
control=control.T


rz=[]
for i in range(len(tename)):
    if len(control[i])>0:
        q75, q25 = np.nanpercentile(control[i], [75 ,25])
        iqr = q75 - q25
        rz.append((score_TE[i]-np.nanmedian(control[i]))/(iqr/1.3489))#robust zscore
    else:
        rz.append(None)



f = open('denovo_'+str(motif)+'_rz_score_enrichment.txt', 'w')
for x in rz:
    f.write(str(x) + "\n")
f.close()
