#!/usr/bin/env python
# coding: utf-8



import csv
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-tf','--tf',type=str,default='N',help='tfname')
parser.add_argument('-pm','--pm',type=str,default='N',help='p:peak,m:,motif')
args=parser.parse_args()
tf=args.tf
pm=args.pm

with open('data/data_3_bigdata_mm9_onlyTE.bed') as f:
    reader = csv.reader(f, delimiter='\t')
    TE_peak = [row for row in reader]


f = open('data/te_tf_mat_TElabel_cortex.txt')
data3 = f.read()
f.close()
tename= data3.split('\n')



tename=list(tename)


TEmat=[[0 for i in range(3)] for j in range(len(tename))]
for i in range(len(TE_peak)):
    TEmat[tename.index(TE_peak[i][3])][0]+=float(((int(TE_peak[i][2])-int(TE_peak[i][1]))))


for i in range(len(TEmat)):
    TEmat[i][2]=float(TEmat[i][0])


leng=[]
for i in range(len(TEmat)):
    leng.append(TEmat[i][2]/(2654895218))

if pm=='p':
    with open('data/TE_contain_'+str(tf)+'_peak_neural_merge.bed') as f:
        reader = csv.reader(f, delimiter='\t')
        tetf = [row for row in reader]
    a=[0]*len(tename)
    a_c=[0]*len(tename)
    for j in range(len(tename)):
        count=[]
        for i in range(len(tetf)):
            if tetf[i][3]==tename[j]:
                a[j]+=min(int(tetf[i][2]),int(tetf[i][-1]))-max(int(tetf[i][1]),int(tetf[i][-2]))
                count.append([tetf[i][0],tetf[i][1],tetf[i][2]])
        a_c[j]+=len(list(map(list, set(map(tuple, count)))))

    out=[]
    for i in range(len(a_c)):
        out.append(a_c[i])
    with open('data/TF_chipatlas_'+str(tf)+'_merge.bed') as f: 
        reader = csv.reader(f, delimiter='\t')
        tfcon = [row for row in reader]
    a_2=0
    for i in range(len(tfcon)):
        a_2+=int(tfcon[i][2])-int(tfcon[i][1])
    
    for i in range(len(a)):
        a[i]=float(a[i])/float(a_2)
    import math
    ans=[]
    for i in range(len(a)):
        if leng[i]==0.0 or a[i]==0.0:
            ans.append('None')
        else:
            ans.append(math.log(a[i]/leng[i],2))

    f = open('out/'+str(tf)+'_TE_log2fold_bp_peak.txt', 'w')
    for x in ans:
        f.write(str(x) + "\n")
    f.close()
    f = open('out/'+str(tf)+'_TE_log2fold_bp_peak_outTE.txt', 'w')
    for x in out:
        f.write(str(x) + "\n")
    f.close()
if pm=='m':
    with open('out/'+str(tf)+'_motif_pos_TEsubfamily_ATAC.txt') as f: 
        reader = csv.reader(f, delimiter='\t')
        motifTE = [row for row in reader]
    with open('out/'+str(tf)+'_motif_pos_ATAC_merge.bed') as f: 
        reader = csv.reader(f, delimiter='\t')
        motifATAC = [row for row in reader]

    motifind=[]
    ATACind=[]
    for i in range(len(motifATAC)):
        motifind.append([motifATAC[i][0],motifATAC[i][1],motifATAC[i][2]])
        ATACind.append([motifATAC[i][-3],motifATAC[i][-2],motifATAC[i][-1]])
    

    a=[0]*len(tename)
    for j in range(len(tename)):
        count=[]
        for i in range(len(motifTE)):
            if motifTE[i][-2]==tename[j]:
                ind=motifind.index([motifTE[i][0],motifTE[i][1],motifTE[i][2]])
                count.append([motifATAC[ind][-3],motifATAC[ind][-2],motifATAC[ind][-1]])
        a[j]+=len(list(map(list, set(map(tuple, count)))))

    out=[]
    for i in range(len(a)):
        out.append(a[i])
    a_2=len(list(map(list, set(map(tuple, ATACind)))))
    for i in range(len(a)):
        a[i]=float(a[i])/float(a_2)
    import math
    ans=[]

    for i in range(len(a)):
        if leng[i]==0.0 or a[i]==0.0:
            ans.append('None')
        else:
            ans.append(math.log(a[i]/leng[i],2))    

    f = open('out/'+str(tf)+'_TE_log2fold_count_motif.txt', 'w')
    for x in ans:
        f.write(str(x) + "\n")
    f.close()
    f = open('out/'+str(tf)+'_TE_log2fold_count_motif_outTE.txt', 'w')
    for x in out:
        f.write(str(x) + "\n")
    f.close()
