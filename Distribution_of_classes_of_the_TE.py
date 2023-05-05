#!/usr/bin/env python
# coding: utf-8

# In[1]:


import csv
import numpy as np
import pandas as pd
import seaborn as sns
import itertools
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import scipy.stats as st
from pybedtools import BedTool
import pybedtools

with open('data/data_3_bigdata_mm9_onlyTE.bed') as f:
    reader = csv.reader(f, delimiter='\t')
    TE = [row for row in reader]
                                                     

cl=['LINE','SINE','LTR','DNA','Unknown']
g_reg=[0,0,0,0,0]
j=0
for j in range(len(cl)):
    tereg=[]
    for i in range(len(TE)):
        if TE[i][-1].split('/')[0]==str(cl[j]):
            tereg.append([TE[i][0],int(TE[i][1]),int(TE[i][2])])
    with open('out/TE_merge_'+str(cl[j])+'.bed','w') as file:
        writer = csv.writer(file,delimiter='\t')
        writer.writerows(tereg)
a = pybedtools.example_bedtool('out/TE_merge_'+str(cl[j])+'.bed')
b=a.sort()
c = b.merge()
reg=0
for i in range(len(c)):
    reg+=c[i].end-c[i].start
g_reg[j]=reg



## ChIP-seq peak

with open('out/TF_chipatlas_Neurod2_merge_TE.bed') as f:
    reader = csv.reader(f, delimiter='\t')
    TEATAC = [row for row in reader]

refcl=['LINE','SINE','LTR','DNA','the_others']
cl=[[],[],[],[],[]]
ans=[0,0,0,0,0]
for i in range(len(TEATAC)):
    if refcl.count(TEATAC[i][-1].split('/')[0])>0:
        cl[refcl.index(TEATAC[i][-1].split('/')[0])].append([TEATAC[i][0],max(int(TEATAC[i][1]),int(TEATAC[i][4])),min(int(TEATAC[i][2]),int(TEATAC[i][5]))])
    else:
        cl[4].append([TEATAC[i][0],max(int(TEATAC[i][1]),int(TEATAC[i][4])),min(int(TEATAC[i][2]),int(TEATAC[i][5]))])

for j in range(len(refcl)):
    with open('out/TF_chipatlas_Neurod2_merge_'+str(refcl[j])+'.bed','w') as file:
        writer = csv.writer(file,delimiter='\t')
        writer.writerows(cl[j])
    a = pybedtools.example_bedtool('out/TF_chipatlas_Neurod2_merge_'+str(refcl[j])+'.bed')
    b=a.sort()
    c = b.merge()
    reg=0
    for i in range(len(c)):
        reg+=c[i].end-c[i].start
    ans[j]=reg
de1reg=ans





with open('out/TF_chipatlas_Lhx2_merge_TE.bed') as f:
    reader = csv.reader(f, delimiter='\t')
    TEATAC = [row for row in reader]

refcl=['LINE','SINE','LTR','DNA','the_others']
cl=[[],[],[],[],[]]
ans=[0,0,0,0,0]
for i in range(len(TEATAC)):
    if refcl.count(TEATAC[i][-1].split('/')[0])>0:
        cl[refcl.index(TEATAC[i][-1].split('/')[0])].append([TEATAC[i][0],max(int(TEATAC[i][1]),int(TEATAC[i][4])),min(int(TEATAC[i][2]),int(TEATAC[i][5]))])
    else:
        cl[4].append([TEATAC[i][0],max(int(TEATAC[i][1]),int(TEATAC[i][4])),min(int(TEATAC[i][2]),int(TEATAC[i][5]))])

for j in range(len(refcl)):
    with open('out/TF_chipatlas_Lhx2_merge_'+str(refcl[j])+'.bed','w') as file:
        writer = csv.writer(file,delimiter='\t')
        writer.writerows(cl[j])
    a = pybedtools.example_bedtool('out/TF_chipatlas_Lhx2_merge_'+str(refcl[j])+'.bed')
    b=a.sort()
    c = b.merge()
    reg=0
    for i in range(len(c)):
        reg+=c[i].end-c[i].start
    ans[j]=reg
de3reg=ans


plt.style.use('default')   
plt.rcParams['savefig.dpi']=400
N = 3
A = np.array([g_reg[0]/sum(g_reg),de1reg[0]/sum(de1reg),de3reg[0]/sum(de3reg)])
B = np.array([g_reg[1]/sum(g_reg),de1reg[1]/sum(de1reg),de3reg[1]/sum(de3reg)])
C = np.array([g_reg[2]/sum(g_reg),de1reg[2]/sum(de1reg),de3reg[2]/sum(de3reg)])
D = np.array([g_reg[3]/sum(g_reg),de1reg[3]/sum(de1reg),de3reg[3]/sum(de3reg)])
E = np.array([g_reg[4]/sum(g_reg),de1reg[4]/sum(de1reg),de3reg[4]/sum(de3reg)])
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence

ind_p = ind + width/2
ind_m = ind - width/2
ind_line = np.sort(np.concatenate([ind_p,ind_m]))
A_line = (np.insert(A, np.arange(3), A))
B_line = (np.insert(B, np.arange(3), B)) + A_line
C_line = (np.insert(C, np.arange(3), C)) + B_line
D_line = (np.insert(D, np.arange(3), D)) + C_line
E_line = (np.insert(E, np.arange(3), E)) + D_line

#plt.figure(figsize=(5,5),dpi=400) 
fig, ax = plt.subplots(dpi=400)

p1 = ax.bar(ind, A, width,zorder=2,hatch='||')
p2 = ax.bar(ind, B, width,bottom=A,zorder=2,hatch='..')
p3 = ax.bar(ind, C, width,bottom=A+B,zorder=2,hatch='x')
p4 = ax.bar(ind, D, width,bottom=A+B+C,zorder=2,hatch='//')
p5 = ax.bar(ind, E, width,bottom=A+B+C+D,zorder=2)


ax.plot(ind_line,A_line,'--k',zorder=1)
ax.plot(ind_line,B_line,'--k',zorder=1)
ax.plot(ind_line,C_line,'--k',zorder=1)
ax.plot(ind_line,D_line,'--k',zorder=1)
ax.plot(ind_line,E_line,'--k',zorder=1)

plt.ylabel('Relative proportion of each TE class')
#plt.title('Scores by group')
#plt.xticks(ind, ('G1', 'G2', 'G3', 'G4'))
plt.xticks([0,1,2], ['genome','Neurod2','Lhx2'])
plt.legend((p1[0], p2[0],p3[0],p4[0],p5[0]), ("LINE", "SINE","LTR","DNA","the others"),
           bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=.1, fontsize=14)
#plt.savefig("ChIP_TE_rate.png",bbox_inches = 'tight', pad_inches = 0)

plt.show()


## ATAC-seq peak

with open('data/merged_cortex_500bp_TE.bed') as f:
    reader = csv.reader(f, delimiter='\t')
    TEATAC = [row for row in reader]





tereg=[[],[],[],[],[]]
for i in range(len(TEATAC)):
    if TEATAC[i][-1].split('/')[0]=='LINE':
        tereg[0].append(TEATAC[i][3])
    elif TEATAC[i][-1].split('/')[0]=='SINE':
        tereg[1].append(TEATAC[i][3])
    elif TEATAC[i][-1].split('/')[0]=='LTR':
        tereg[2].append(TEATAC[i][3])
    elif TEATAC[i][-1].split('/')[0]=='DNA':
        tereg[3].append(TEATAC[i][3])
    else:
        tereg[4].append(TEATAC[i][3])                                             


leng=[0,0,0,0,0]
for i in range(len(tereg)):
    leng[i]=len(list(set(tereg[i])))


tereg=leng





with open('out/denovo_1_motif_pos_TEsubfamily_ATAC.bed') as f:
    reader = csv.reader(f, delimiter='\t')
    motifTE = [row for row in reader]

refcl=['LINE','SINE','LTR','DNA','the others']
cl=[[],[],[],[],[]]
for i in range(len(motifTE)):
    if refcl.count(motifTE[i][-1].split('/')[0])>0:
        cl[refcl.index(motifTE[i][-1].split('/')[0])].append('_'.join(motifTE[i][5:8]))
    else:
        cl[4].append('_'.join(motifTE[i][5:8]))
leng=[]
for i in range(len(cl)):
    leng.append(len(list(set(cl[i]))))

de1reg=leng


with open('out/denovo_2_motif_pos_TEsubfamily_ATAC.bed') as f:
    reader = csv.reader(f, delimiter='\t')
    motifTE = [row for row in reader]

refcl=['LINE','SINE','LTR','DNA','the others']
cl=[[],[],[],[],[]]
for i in range(len(motifTE)):
    if refcl.count(motifTE[i][-1].split('/')[0])>0:
        cl[refcl.index(motifTE[i][-1].split('/')[0])].append('_'.join(motifTE[i][5:8]))
    else:
        cl[4].append('_'.join(motifTE[i][5:8]))
leng=[]
for i in range(len(cl)):
    leng.append(len(list(set(cl[i]))))

de2reg=leng


with open('out/denovo_3_motif_pos_TEsubfamily_ATAC.bed') as f:
    reader = csv.reader(f, delimiter='\t')
    motifTE = [row for row in reader]

refcl=['LINE','SINE','LTR','DNA','the others']
cl=[[],[],[],[],[]]
for i in range(len(motifTE)):
    if refcl.count(motifTE[i][-1].split('/')[0])>0:
        cl[refcl.index(motifTE[i][-1].split('/')[0])].append('_'.join(motifTE[i][5:8]))
    else:
        cl[4].append('_'.join(motifTE[i][5:8]))
leng=[]
for i in range(len(cl)):
    leng.append(len(list(set(cl[i]))))

de3reg=leng






plt.style.use('default') 
plt.rcParams['savefig.dpi']=400
N = 5
A = np.array([g_reg[0]/sum(g_reg),tereg[0]/sum(tereg),de1reg[0]/sum(de1reg),de2reg[0]/sum(de2reg),de3reg[0]/sum(de3reg)])
B = np.array([g_reg[1]/sum(g_reg),tereg[1]/sum(tereg),de1reg[1]/sum(de1reg),de2reg[1]/sum(de2reg),de3reg[1]/sum(de3reg)])
C = np.array([g_reg[2]/sum(g_reg),tereg[2]/sum(tereg),de1reg[2]/sum(de1reg),de2reg[2]/sum(de2reg),de3reg[2]/sum(de3reg)])
D = np.array([g_reg[3]/sum(g_reg),tereg[3]/sum(tereg),de1reg[3]/sum(de1reg),de2reg[3]/sum(de2reg),de3reg[3]/sum(de3reg)])
E = np.array([g_reg[4]/sum(g_reg),tereg[4]/sum(tereg),de1reg[4]/sum(de1reg),de2reg[4]/sum(de2reg),de3reg[4]/sum(de3reg)])
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence

ind_p = ind + width/2
ind_m = ind - width/2
ind_line = np.sort(np.concatenate([ind_p,ind_m]))
A_line = (np.insert(A, np.arange(5), A))
B_line = (np.insert(B, np.arange(5), B)) + A_line
C_line = (np.insert(C, np.arange(5), C)) + B_line
D_line = (np.insert(D, np.arange(5), D)) + C_line
E_line = (np.insert(E, np.arange(5), E)) + D_line

#plt.figure(figsize=(5,5),dpi=400) 
fig, ax = plt.subplots(dpi=400)

p1 = ax.bar(ind, A, width,zorder=2,hatch='||')
p2 = ax.bar(ind, B, width,bottom=A,zorder=2,hatch='..')
p3 = ax.bar(ind, C, width,bottom=A+B,zorder=2,hatch='x')
p4 = ax.bar(ind, D, width,bottom=A+B+C,zorder=2,hatch='//')
p5 = ax.bar(ind, E, width,bottom=A+B+C+D,zorder=2)

ax.plot(ind_line,A_line,'--k',zorder=1)
ax.plot(ind_line,B_line,'--k',zorder=1)
ax.plot(ind_line,C_line,'--k',zorder=1)
ax.plot(ind_line,D_line,'--k',zorder=1)
ax.plot(ind_line,E_line,'--k',zorder=1)

plt.ylabel('Relative proportion of each TE class')
#plt.title('Scores by group')
#plt.xticks(ind, ('G1', 'G2', 'G3', 'G4'))
plt.xticks([0,1, 2,3,4], ['genome','ATAC peak \n (n=63059)','ATAC peak \n with $\it{de}$ $\it{novo}$ 1 \n (n=4744)','ATAC peak \n with $\it{de}$ $\it{novo}$ 2 \n (n=4356)','ATAC peak \n with $\it{de}$ $\it{novo}$ 3 \n (n=930)'],rotation=30)
plt.legend((p1[0], p2[0],p3[0],p4[0],p5[0]), ("LINE", "SINE","LTR","DNA","the others"),
           bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=.1, fontsize=14)
#plt.savefig("ATAC_TE_rate.png",bbox_inches = 'tight', pad_inches = 0)

plt.show()


