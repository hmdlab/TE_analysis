#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-te','--te',type=str,default='None',help='TE name')
parser.add_argument('-tf','--tf',type=str,default='None',help='TF ChIP name')
args=parser.parse_args()
te=args.te 
tf=args.tf 


# In[18]:


import numpy as np
import csv
import gzip
import re


with gzip.open('data/mm9.fa.align.gz','rt') as f:
    data3 = f.read()  
    f.close()
    align= data3.split('\n')

ch=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
for i in range(len(align)):
    align[i]=align[i].split(" ")

## mapping from genome to reference TE

    
with open('data/mm9.fa.align_map_onlyTE.bed') as f: 
    reader = csv.reader(f, delimiter='\t')
    ref = [row for row in reader]


reflabel=[]
telabel=[]
for i in range(len(ref)):
    reflabel.append([ref[i][0],ref[i][1],ref[i][2],ref[i][3]])
    telabel.append(ref[i][3])


with open('data/'+str(tf)+'_peak_mm9.fa.align_map_only_TE.bed') as f:
    reader = csv.reader(f, delimiter='\t')
    tetf = [row for row in reader]
tetf_c=[]
for i in range(len(tetf)):
    tetf_c.append([tetf[i][0],tetf[i][1],tetf[i][2],tetf[i][3],tetf[i][4],tetf[i][5],tetf[i][6],tetf[i][7],tetf[i][-3],tetf[i][-2],tetf[i][-1],str(tf)])
tetf=tetf_c
ok=[]
for i in range(len(tetf)):
    if tetf[i][3]==te and tetf[i][-1]==tf:
        ok.append(tetf[i])
tetf=ok

index=reflabel.index([tetf[0][0],tetf[0][1],tetf[0][2],te])
if len(ref[index][5].split('('))>1:
    rang=int(re.split('[()]',ref[index][5])[1])+int(ref[index][6])
elif len(ref[index][7].split('('))>1:
    rang=int(re.split('[()]',ref[index][7])[1])+int(ref[index][6])

score_range=[0 for i in range(rang)]#count
control=[0 for i in range(rang)]#control
for i in range(len(tetf)):
    index=reflabel.index([tetf[i][0],tetf[i][1],tetf[i][2],te])
    #print(ref[index])
    if float(ref[index][-2])==0.0 and float(ref[index][-3])==0.0:
        if len(ref[index][5].split('('))>1:#reverse
            for j in range(int(ref[index][-5])-1-max(int(tetf[i][1]),int(tetf[i][9]))+int(tetf[i][1]),int(ref[index][-4])-2-min(int(tetf[i][2]),int(tetf[i][10]))+int(tetf[i][2]),-1):#
                score_range[j]+=1
        elif len(ref[index][7].split('('))>1:#forward
            for j in range(int(ref[index][-6])-1+max(int(tetf[i][1]),int(tetf[i][9]))-int(tetf[i][1]),int(ref[index][-5])+min(int(tetf[i][2]),int(tetf[i][10]))-int(tetf[i][2])):
                score_range[j]+=1
    else:#insert or deletion
        if len(ref[index][5].split('('))>1:
            k=int(ref[index][-1])+4
            length=len([a for a in align[k-4] if a != ''][5])
            bp_del=[]
            while True:
                list_a = [a for a in align[k] if a != '']
                list_g = [a for a in align[k-2] if a != '']
                if list_a[0]=='C':
                    if len(list_g[0].split(':'))>1:
                        now=len(list_g[1])
                        sb=int(list_g[0].split(':')[1].split('-')[0]+'0'*(length-len(list_g[0].split(':')[1].split('-')[0])))+int(list_g[1])
                    else:
                        sb=int(list_g[1])
                    #print(list_g)
                    #print(list_a)
                    enter=[n for n, v in enumerate(list_a[3]) if v == '-']
                    leave=[n for n, v in enumerate(list_g[2]) if v == '-']
                    for e in range(len(enter)):
                        bp_del.append(sb+enter[e]-np.count_nonzero(np.array(leave)<enter[e]))
                else:
                    if len(list_g[0].split(':'))>1:
                        now=len(list_g[1])
                        sb=int(list_g[0].split(':')[1].split('-')[0]+'0'*(length-len(list_g[0].split(':')[1].split('-')[0])))+int(list_g[1])
                    else:
                        sb=int(list_g[1])
                    #print(list_g)
                    #print(list_a)
                    enter=[n for n, v in enumerate(list_a[2]) if v == '-']
                    leave=[n for n, v in enumerate(list_g[2]) if v == '-']
                    for e in range(len(enter)):
                        bp_del.append(sb+enter[e]-np.count_nonzero(np.array(leave)<enter[e]))
                if [a for a in align[k] if a != ''][-1]==ref[index][7]:
                    break
                k+=4
            b_del=[]
            out_b=[]
            out_o=[]
            for d in range(len(bp_del)):
                bp_del_in=bp_del[d]-int(tetf[i][9])
                if bp_del_in>=0 and bp_del_in<int(tetf[i][10])-int(tetf[i][9])+1:
                    b_del.append(bp_del_in)
                elif bp_del_in<0:
                    out_b.append(bp_del_in)
                elif bp_del_in>=int(tetf[i][10])-int(tetf[i][9])+1:
                    out_o.append(bp_del_in)
            if float(ref[index][-3])!=0.0:#deletion
                k=int(ref[index][-1])+2
                bp_in=[]
                gene_in=[]
                while True:
                    gene_in_b=[]
                    list_g = [a for a in align[k] if a != '']
                    list_a = [a for a in align[k+2] if a != '']
                    #print(list_a)
                    #print(list_g)
                    if list_a[0]=='C':
                        now=len(list_g[1])
                        sb=int(list_a[2])
                        enter=[n for n, v in enumerate(list_a[3]) if v == '-']
                        leave=[n for n, v in enumerate(list_g[2]) if v == '-']
                        if len(list_g[0].split(':'))>1:
                            gene=int(list_g[0].split(':')[1].split('-')[0]+'0'*(length-len(list_g[0].split(':')[1].split('-')[0])))+int(list_g[1])
                        else:
                            gene=int(list_g[1])
                    else:
                        now=len(list_g[1])
                        sb=int(list_a[1])
                        enter=[n for n, v in enumerate(list_a[2]) if v == '-']
                        leave=[n for n, v in enumerate(list_g[2]) if v == '-']
                        if len(list_g[0].split(':'))>1:
                            gene=int(list_g[0].split(':')[1].split('-')[0]+'0'*(length-len(list_g[0].split(':')[1].split('-')[0])))+int(list_g[1])
                        else:
                            gene=int(list_g[1])
                    gene_in_b=[n+gene for n, v in enumerate(list_g[2]) if v == '-']
                    for g in range(len(gene_in_b)):
                        gene_in_b[g]=gene_in_b[g]-g-1
                    gene_in+=gene_in_b
                    for e in range(len(leave)):
                        bp_in.append(sb-leave[e]+np.count_nonzero(np.array(enter)<leave[e])-1)
                    if [a for a in align[k+2] if a != ''][-1]==ref[index][7]:
                        break
                    k+=4
                in_b=[]
                in_o=[]
                for m in range(len(gene_in)):
                    if gene_in[m]<int(tetf[i][9]):
                        in_b.append(m)
                    elif gene_in[m]>int(tetf[i][10])-1:
                        in_o.append(m)
                #print('in_b',in_b,'in_o',in_o,'out_b',out_b,'out_o',out_o)
                for j in range(int(ref[index][-5])-1-max(int(tetf[i][1]),int(tetf[i][9]))+int(tetf[i][1])-len(in_b)+len(out_b),int(ref[index][-4])-2-min(int(tetf[i][2]),int(tetf[i][10]))+int(tetf[i][2])+len(in_o)-len(out_o),-1):#
                    if bp_in.count(j)==0:
                        score_range[j]+=1
            else:
                for j in range(int(ref[index][-5])-1-max(int(tetf[i][1]),int(tetf[i][9]))+int(tetf[i][1])+len(out_b),int(ref[index][-4])-2-min(int(tetf[i][2]),int(tetf[i][10]))+int(tetf[i][2])-len(out_o),-1):#
                    score_range[j]+=1
        elif len(ref[index][7].split('('))>1:#forward
            k=int(ref[index][-1])+4
            length=len([a for a in align[k-4] if a != ''][5])
            bp_del=[]
            while True:
                list_a = [a for a in align[k] if a != '']
                list_g = [a for a in align[k-2] if a != '']
                if list_a[0]=='C':
                    if len(list_g[0].split(':'))>1:
                        now=len(list_g[1])
                        sb=int(list_g[0].split(':')[1].split('-')[0]+'0'*(length-len(list_g[0].split(':')[1].split('-')[0])))+int(list_g[1])
                    else:
                        sb=int(list_g[1])
                    #print(list_g)
                    #print(list_a)
                    enter=[n for n, v in enumerate(list_a[3]) if v == '-']
                    leave=[n for n, v in enumerate(list_g[2]) if v == '-']
                    for e in range(len(enter)):
                        bp_del.append(sb+enter[e]-np.count_nonzero(np.array(leave)<enter[e]))
                else:
                    if len(list_g[0].split(':'))>1:
                        now=len(list_g[1])
                        sb=int(list_g[0].split(':')[1].split('-')[0]+'0'*(length-len(list_g[0].split(':')[1].split('-')[0])))+int(list_g[1])
                    else:
                        sb=int(list_g[1])
                    #print(list_g)
                    #print(list_a)
                    enter=[n for n, v in enumerate(list_a[2]) if v == '-']
                    leave=[n for n, v in enumerate(list_g[2]) if v == '-']
                    for e in range(len(enter)):
                        bp_del.append(sb+enter[e]-np.count_nonzero(np.array(leave)<enter[e]))
                if int([a for a in align[k] if a != ''][-1])==int(ref[index][6]):
                    break
                #print('ref:',ref[index][6],int([a for a in align[k] if a != ''][-1]))
                k+=4
            b_del=[]
            out_b=[]
            out_o=[]
            for d in range(len(bp_del)):
                bp_del_in=bp_del[d]-int(tetf[i][9])
                if bp_del_in>=0 and bp_del_in<int(tetf[i][10])-int(tetf[i][9])+1:
                    b_del.append(bp_del_in)
                elif bp_del_in<0:
                    out_b.append(bp_del_in)
                elif bp_del_in>=int(tetf[i][10])-int(tetf[i][9])+1:
                    out_o.append(bp_del_in)
            if float(ref[index][-3])!=0.0:#deletion
                k=int(ref[index][-1])+2
                bp_in=[]
                gene_in=[]
                while True:
                    gene_in_b=[]
                    list_g = [a for a in align[k] if a != '']
                    list_a = [a for a in align[k+2] if a != '']
                    #print(list_a)
                    #print(list_g)
                    if list_a[0]=='C':
                        sb=int(list_a[2])
                        enter=[n for n, v in enumerate(list_a[3]) if v == '-']
                        leave=[n for n, v in enumerate(list_g[2]) if v == '-']
                        if len(list_g[0].split(':'))>1:
                            now=len(list_g[1])
                            gene=int(list_g[0].split(':')[1].split('-')[0]+'0'*(length-len(list_g[0].split(':')[1].split('-')[0])))+int(list_g[1])
                        else:
                            gene=int(list_g[1])
                    else:
                        sb=int(list_a[1])
                        enter=[n for n, v in enumerate(list_a[2]) if v == '-']
                        leave=[n for n, v in enumerate(list_g[2]) if v == '-']
                        if len(list_g[0].split(':'))>1:
                            now=len(list_g[1])
                            gene=int(list_g[0].split(':')[1].split('-')[0]+'0'*(length-len(list_g[0].split(':')[1].split('-')[0])))+int(list_g[1])
                        else:
                            gene=int(list_g[1])
                    gene_in_b=[n+gene for n, v in enumerate(list_g[2]) if v == '-']
                    for g in range(len(gene_in_b)):
                        gene_in_b[g]=gene_in_b[g]-g-1
                    gene_in+=gene_in_b
                    for e in range(len(leave)):
                        bp_in.append(sb+leave[e]-np.count_nonzero(np.array(enter)<leave[e])-1)
                    if [a for a in align[k+2] if a != ''][-1]==ref[index][6]:
                        break
                    k+=4
                in_b=[]
                in_o=[]
                for m in range(len(gene_in)):
                    if gene_in[m]<int(tetf[i][9]):
                        in_b.append(m)
                    elif gene_in[m]>int(tetf[i][10])-1:
                        in_o.append(m)
                for j in range(int(ref[index][-6])-1+max(int(tetf[i][1]),int(tetf[i][9]))-int(tetf[i][1])+len(in_b)-len(out_b),int(ref[index][-5])+min(int(tetf[i][2]),int(tetf[i][10]))-int(tetf[i][2])-len(in_o)+len(out_o)):
                    if bp_in.count(j)==0:
                        score_range[j]+=1  
            else:
                for j in range(int(ref[index][-6])-1+max(int(tetf[i][1]),int(tetf[i][9]))-int(tetf[i][1])-len(out_b),int(ref[index][-5])+min(int(tetf[i][2]),int(tetf[i][10]))-int(tetf[i][2])+len(out_o)):
                    score_range[j]+=1

f = open('out/'+str(tf)+'_'+str(te)+'.Freq_refTE_count.txt', 'w')
for x in score_range:
    f.write(str(x) + "\n")
f.close()

for i in range(len(telabel)):
    if telabel[i]==te:
        if float(ref[i][-3])==0.0:
            if len(ref[i][5].split('('))>1:#reverse
                for j in range(int(ref[i][-5])-1,int(ref[i][-4])-2,-1):
                    control[j]+=1
            elif len(ref[i][7].split('('))>1:#forward
                for j in range(int(ref[i][-6])-1,int(ref[i][-5])):
                    control[j]+=1
        if float(ref[i][-3])!=0.0:#deletion
            if len(ref[i][5].split('('))>1:
                k=int(ref[i][-1])+2
                bp_in=[]
                while True:
                    if [a for a in align[k] if a != ''][0]=='Matrix':
                        break
                    list_g = [a for a in align[k] if a != '']
                    list_a = [a for a in align[k+2] if a != '']
                    if list_a[0]=='C':
                        sb=int(list_a[2])
                        enter=[n for n, v in enumerate(list_a[3]) if v == '-']
                        leave=[n for n, v in enumerate(list_g[2]) if v == '-']
                    else:
                        sb=int(list_a[1])
                        enter=[n for n, v in enumerate(list_a[2]) if v == '-']
                        leave=[n for n, v in enumerate(list_g[2]) if v == '-']
                    for e in range(len(leave)):
                        bp_in.append(sb-leave[e]+np.count_nonzero(np.array(enter)<leave[e])-1)
                    #print(int([a for a in align[k+2] if a != ''][-1]),int(ref[i][7]))
                    if int([a for a in align[k+2] if a != ''][-1])==int(ref[i][7]):
                        break
                    k+=4
                for j in range(int(ref[i][-5])-1,int(ref[i][-4])-2,-1):
                    if bp_in.count(j)==0:
                        control[j]+=1
            elif len(ref[i][7].split('('))>1:#forward
                k=int(ref[i][-1])+2
                bp_in=[]
                while True:
                    if [a for a in align[k] if a != ''][0]=='Matrix':
                        break
                    list_g = [a for a in align[k] if a != '']
                    list_a = [a for a in align[k+2] if a != '']
                    if list_a[0]=='C':
                        sb=int(list_a[2])
                        enter=[n for n, v in enumerate(list_a[3]) if v == '-']
                        leave=[n for n, v in enumerate(list_g[2]) if v == '-']
                    else:
                        sb=int(list_a[1])
                        enter=[n for n, v in enumerate(list_a[2]) if v == '-']
                        leave=[n for n, v in enumerate(list_g[2]) if v == '-']
                    for e in range(len(leave)):
                        bp_in.append(sb+leave[e]-np.count_nonzero(np.array(enter)<leave[e])-1)
                    #print(int([a for a in align[k+2] if a != ''][-1]),int(ref[i][6]))
                    if ref[i][5]==ref[i][6]:
                        break
                    if int([a for a in align[k+2] if a != ''][-1])==int(ref[i][6]):
                        break
                    k+=4          
                for j in range(int(ref[i][-6])-1,int(ref[i][-5])):
                    if bp_in.count(j)==0:
                        control[j]+=1

f = open('out/'+str(tf)+'_'+str(te)+'.Freq_refTE_control.txt', 'w')
for x in control:
    f.write(str(x) + "\n")
f.close()
