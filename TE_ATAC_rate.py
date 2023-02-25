#!/usr/bin/env python
# coding: utf-8

import csv
with open('out/cell_type_label.txt') as f:
    reader = csv.reader(f, delimiter=' ')
    celltype = [row for row in reader]
copy=[]
for i in range(len(celltype)):
    copy.append(celltype[i][1])
celltype=copy

from scipy.sparse import csr_matrix, csc_matrix, coo_matrix, lil_matrix
import scipy.io
mat = scipy.io.mmread("out/counts_filtered.mtx")
matrix= mat.tocsr()
matrix=matrix.toarray()


f = open('out/celllabel.txt')
data3 = f.read() 
f.close()
cellid= data3.split('\n')

with open('out/merged_cortex_500bp.txt') as f:
    reader = csv.reader(f, delimiter='\t')
    peak = [row for row in reader]

peakid=[]
for i  in range(len(peak)):
    peakid.append("_".join(peak[i][:3])



df = pd.DataFrame(matrix.T,columns=list(peakid[:-1]),index=list(cellid[:-1]))


df_b=df.where(df > 0, 1)


df_b = pd.DataFrame(np.where(df > 0, 1, df),
                           index=df.index, columns=df.columns)

sum_mat=df_b.sum(axis=1)



with open('out/merged_cortex_500bp_TE.bed') as f: 
    reader = csv.reader(f, delimiter='\t')
    TE = [row for row in reader]


TEr=[]
for i in range(len(TE)):
    TEr.append("_".join(TE[i][:3])
TEr=list(set(TEr))



TE_a=df_b.T.filter(items=TEr, axis='index')


TE_a=TE_a.T



sum_TEmat=TE_a.sum(axis=1)

label=['Ex._neurons_CPN','Ex._neurons_SCPN','Ex._neurons_CThPN','SOM+_Interneurons','Inhibitory_neurons',
         'Oligodendrocytes','Microglia','Astrocytes']
te=[]
for j in range(len(label)):
    celltype_ind=[]
    for i in range(len(celltype)):
        if celltype[i]==label[j]:
            celltype_ind.append(sum_mat.index[i])
    sum_mat_cell=sum_mat.filter(items=celltype_ind, axis='index')
    sum_TEmat_cell=sum_TEmat.filter(items=celltype_ind, axis='index')
    ans=[]
    for i in range(len(celltype_ind)):
        ans.append((sum_TEmat_cell/sum_mat_cell)[i])
    te.append(ans)


tec=te
for i in range(len(te)):
    for j in range(len(te[i])):
        te[i][j]=float(tec[i][j])*100


plt.rcParams["font.size"] = 18


from matplotlib import pyplot as plt
import seaborn as sns
name=label
plt.style.use('default')
sns.set()
sns.set_style('whitegrid')
sns.set_palette('gray')

heights = []
std = []
for i in range(len(te)):
    heights.append(np.mean(te[i])+0.01)
    std.append(np.std(te[i]))
width = 0.8 # the width of the bars
bars = np.arange(len(heights))
    
fig = plt.figure(figsize=(8,8),dpi=400)
ax = fig.add_subplot(1, 1, 1)

ax.violinplot(te,showmeans=True)
ax.set_xticks([1, 2, 3,4,5,6,7,8])
ax.set_xticklabels(name,rotation=30,fontsize=14)

ax.set_xlabel('cell type')
ax.set_ylabel('TE-derived ATAC-seq peak count /\nATAC-seq peak count\n[%]')
plt.show()
