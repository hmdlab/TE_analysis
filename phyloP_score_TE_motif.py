#!/usr/bin/env python
# coding: utf-8


import csv
import numpy as np


import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-te','--te',type=str,default='None',help='TEname')
parser.add_argument('-tf','--tf',type=str,default='None',help='TFname')
parser.add_argument('-cont','--cont',type=str,default='None',help='control')
parser.add_argument('-pre','--pre',type=str,default='None',help='prepare Y:yes')
args=parser.parse_args()
te=args.te 
tf=args.tf 
cont=args.cont
pre=args.pre

if pre=='Y':
    with open('data/data_3_bigdata_mm9_onlyTE.bed') as f:
        reader = csv.reader(f, delimiter='\t')
        TE_peak = [row for row in reader]
    
    ans=[]
    for i in range(len(TE_peak)):
        if TE_peak[i][-2]==te:
            ans.append(TE_peak[i])
            
    with open('data/'+str(te)+'_region.bed','w') as file:
        writer = csv.writer(file,delimiter='\t')
        writer.writerows(ans)   
    
    con=[]
    for i in range(1000):
        con.append(['chr1',100000,100006])
        
    with open('data/control_1000_7bp.bed','w') as file:
        writer = csv.writer(file,delimiter='\t')
        writer.writerows(con)   
        
else:
    with open('data/'+str(tf)+'_motif_pos_TEsubfamily_ATAC.txt') as f: 
        reader = csv.reader(f, delimiter='\t')
        TETF_peak = [row for row in reader]


    tetf=[]
    tetf_c=[]
    for i in range(len(TETF_peak)):
        if tetf_c.count([TETF_peak[i][0],TETF_peak[i][1],TETF_peak[i][2]])==0:
            tetf.append(TETF_peak[i])
            tetf_c.append([TETF_peak[i][0],TETF_peak[i][1],TETF_peak[i][2]])
    TETF_peak=tetf



    tetf=[]
    for k in range(len(tf)):
        for i in range(len(TETF_peak)):
            if te=='None':
                tetf.append(TETF_peak[i])
            else:
                if TETF_peak[i][8]==te:
                    tetf.append(TETF_peak[i])



    ch=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
    name=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']


    import gzip
    for i in range(len(ch)):
        with gzip.open('data/'+str(ch[i])+'.phyloP30way.placental.wigFix.gz','rt') as f:
            data3 = f.read()  
            f.close()
            ch[i]= data3.split('\n')


    # In[ ]:


    start_bp=[[0 for i in range(0)] for j in range(len(name))]
    num=[[0 for i in range(0)] for j in range(len(name))]
    for k in range(len(name)):
        for i in range(len(ch[k])):
            if len(ch[k][i].split(" "))>1:
                start_chr=ch[k][i].split(" ")
                start_number=int(start_chr[2].split("=")[1])
                start_bp[k].append(start_number)
                num[k].append(i)

    for i in range(len(name)):
        num[i].append(len(ch[i]))
        start_bp[i].append(start_bp[i][-1]+len(ch[i]))

    score=[]
    score_bp=[]
    if cont=='True':
        with open('out/'+str(te)+'ATAC_intersect_control.bed') as f: 
            reader = csv.reader(f, delimiter='\t')
            tetf = [row for row in reader]
        for k in range(len(tetf)):
            ct=[]
            for i in range(int(tetf[k][1]),int(tetf[k][2])):
                ct.append(float(ch[name.index(tetf[k][0])][i]))
            score.append(sum(ct)/len(ct))

        TEscore=[]
        for i in range(len(score)):
            TEscore.append([tetf[i][0],tetf[i][1],tetf[i][2],tetf[i][-2],score[i]])
        with open('out/'+str(te)+'.phyloPsocre_control.bed','w') as file:
            writer = csv.writer(file,delimiter='\t')
            writer.writerows(score)
    else:
        for k in range(len(tetf)):
            ct=[]
            if tetf[i][-2]==te:
                for i in range(int(tetf[k][1]),int(tetf[k][2])):
                    ct.append(float(ch[name.index(tetf[k][0])][i]))
                score.append(sum(ct)/len(ct))

        TEscore=[]
        for i in range(len(score)):
            TEscore.append([tetf[i][0],tetf[i][1],tetf[i][2],tetf[i][-2],score[i]])
        with open('out/'+str(te)+'.phyloPsocre_'+str(tf)+'.bed','w') as file:
            writer = csv.writer(file,delimiter='\t')
            writer.writerows(score)

    #print(np.average(score))

