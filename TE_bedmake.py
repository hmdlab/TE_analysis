#!/usr/bin/env python
# coding: utf-8

# In[1]:

import csv
nonTE=['Simple_repeat','Low_complexity','Satellite','scRNA','tRNA','snRNA','srpRNA','RNA','rRNA']
ch=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']

with open('data/mm9.fa.out') as f:
    reader = csv.reader(f, delimiter=' ')
    TE = [row for row in reader]

ans=[]
for i in range(3,len(TE)):
    c=[a for a in TE[i] if a != '']
    if ch.count(c[4])>0 and nonTE.count(c[-5])==0:
        ans.append([c[4],c[5],c[6],c[-6],c[-5]])

with open('data/data_3_bigdata_mm9_onlyTE.bed','w') as file:
    writer = csv.writer(file,delimiter='\t')
    writer.writerows(ans)

    
with gzip.open('data/mm9.fa.align.gz','rt') as f:
    data3 = f.read()  
    f.close()
    align= data3.split('\n')

for i in range(len(align)):
    align[i]=align[i].split(" ")

ref=[]
for i in range(len(align)):
    c=[a for a in align[i] if a != '']
    if len(c)>6:
        if ch.count(c[4])>0 and nonTE.count(c[-4])==0:
            ref.append([c[4],c[5],c[6],c[-5],c[-4],c[-3],c[-2],c[-1],c[2],c[3],i])
with open('data/mm9.fa.align_map_onlyTE.bed','w') as file:
    writer = csv.writer(file,delimiter='\t')
    writer.writerows(ref)   
    
