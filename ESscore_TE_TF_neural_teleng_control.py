#!/usr/bin/env python
# coding: utf-8

# In[214]:



import csv
from re import S
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix, lil_matrix
import scipy.io
import collections
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-tf','--tf',type=str,default='N',help='tfname')
parser.add_argument('-pm','--pm',type=str,default='N',help='p:peak,m:,motif')
parser.add_argument('-pre','--pre',type=str,default='N',help='control prepare')
args=parser.parse_args()
tf=args.tf
pm=args.pm
pre=args.pre

with open('data/data_3_bigdata_mm9_onlyTE.bed') as f:
    reader = csv.reader(f, delimiter='\t')
    TE_peak = [row for row in reader]

# In[151]:


tename=[]
for i in range(len(TE_peak)):
    tename.append(TE_peak[i][3])
telabel=list(set(tename))



f = open('data/te_tf_mat_TElabel_cortex.txt')
data3 = f.read()
f.close()
tename= data3.split('\n')


# In[119]:


#tename=np.delete(tename,-1)


tename=list(tename)


TEmat=[[0 for i in range(3)] for j in range(len(tename))]
for i in range(len(TE_peak)):
    TEmat[tename.index(TE_peak[i][3])][0]+=float(((int(TE_peak[i][2])-int(TE_peak[i][1]))))


# In[1872]:


for i in range(len(TEmat)):
    TEmat[i][2]=float(TEmat[i][0])


leng=[]
s3=[]
s4=[]
for i in range(len(TEmat)):
    s3.append(TEmat[i][2])
    s4.append(2654895218)
    leng.append(TEmat[i][2]/(2654895218))

with open('out/background_peaks_kmers_500bp.csv') as f:
    reader = csv.reader(f, delimiter=' ')
    bg = [row for row in reader]

cp=[]
for i in range(1,len(bg)):
    cp.append(bg[i][1:])
bg=cp
bg=np.array(bg,dtype='float')
bg=np.array(bg,dtype='int')

if pm=='m' and pre=='Y':
    f = open('out/merged_cortex_500bp.txt')
    data3 = f.read()
    f.close()
    peakid= data3.split('\n')
    peakid=peakid[:-1]
    peak=[]
    for i in range(len(peakid)):
        peak.append(peakid[i].split('_'))


    with open('out/'+str(tf)+'_motif_pos_ATAC_merge.bed') as f: 
        reader = csv.reader(f, delimiter='\t')
        motifATAC = [row for row in reader]
    
    for i in range(len(bg[0])):  
        cp=[]  
        for  j in range(len(motifATAC)):
            ind_p=peakid.index('_'.join(motifATAC[j][-3:]))
            ind=bg[ind_p,i]-1
            cp.append([peak[ind][0],int(peak[ind][1])+int(motifATAC[j][1])-int(motifATAC[j][-2]),int(peak[ind][1])+int(motifATAC[j][2])-int(motifATAC[j][-2]),motifATAC[j][3],motifATAC[j][4],peak[ind][0],peak[ind][1],peak[ind][2]])
        motifATAC=cp
        with open('out/'+str(tf)+'_motif_pos_ATAC_merge_control_'+str(i+1)+'.bed','w') as file:
            writer = csv.writer(file,delimiter='\t')
            writer.writerows(motifATAC)



elif pm=='m':
    for k in range(len(bg[0])):
        with open('out/'+str(tf)+'_motif_pos_TEsubfamily_ATA_control_'+str(k+1)+'.txt') as f: 
            reader = csv.reader(f, delimiter='\t')
            motifTE = [row for row in reader]
        with open('out/'+str(tf)+'_motif_pos_ATAC_merge_control_'+str(k+1)+'.bed') as f: 
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
        s1=[]
        s2=[]
        for i in range(len(a)):
            out.append(a[i])
        a_2=len(list(map(list, set(map(tuple, ATACind)))))
        for i in range(len(a)):
            s1.append(a[i])
            s2.append(a_2)
            a[i]=float(a[i])/float(a_2)
        import math
        ans=[]
        

        for i in range(len(a)):
            if leng[i]==0.0 or a[i]==0.0:
                ans.append('None')
            else:
                ans.append(math.log(a[i]/leng[i],2))    

        s=[s1,s2,s3,s4]
        with open('out/'+str(tf)+'_TE_log2fold_count_motif_control_'+str(k+1)+'_detail.bed', 'w') as file:
            writer = csv.writer(file,delimiter='\t')
            writer.writerows(s)

        f = open('out/'+str(tf)+'_TE_log2fold_count_motif_control_'+str(k+1)+'.txt', 'w')
        for x in ans:
            f.write(str(x) + "\n")
        f.close()
        f = open('out/'+str(tf)+'_TE_log2fold_count_motif_outTE_control_'+str(k+1)+'.txt', 'w')
        for x in out:
            f.write(str(x) + "\n")
        f.close()
