#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-te','--te',type=str,default='None',help='TEname')
parser.add_argument('-denovo','--denovo',type=str,default='None',help='Y:denovo motif_number')
args=parser.parse_args()
te=args.te 
denovo=args.denovo


# In[18]:


import numpy as np
import csv
import gzip
import re


with gzip.open('data/mm9.fa.align.gz','rt') as f:
    data3 = f.read() 
    f.close()
    align= data3.split('\n')

for i in range(len(align)):
    align[i]=align[i].split(" ")




with open('data/mm9.fa.align_map_onlyTE.bed') as f: 
    reader = csv.reader(f, delimiter='\t')
    ref = [row for row in reader]



reflabel=[]
telabel=[]
for i in range(len(ref)):
    reflabel.append([ref[i][0],ref[i][1],ref[i][2],ref[i][3]])
    telabel.append(ref[i][3])


if denovo!='None':
    with open('out/denovo_'+str(denovo)+'_motif_pos_TEsubfamily_ATAC_align.bed') as f:
        reader = csv.reader(f, delimiter='\t')
        motif_b = [row for row in reader]    

cp=[]
for i in range(len(motif_b)):
    cp.append([motif_b[i][0],motif_b[i][1],motif_b[i][2],motif_b[i][3],motif_b[i][4],motif_b[i][5],motif_b[i][6],motif_b[i][7],motif_b[i][8],motif_b[i][9]])

motif_b=cp

if mask!='N':
    motif_tf=[]
    motif_a=[]
    for i in range(len(motif_b)):
        if len(motif_b[i])>0 and motif_b[i][-2]==te:
            motif_a.append(motif_b[i])
            motif_tf.append([motif_b[i][0],motif_b[i][1],motif_b[i][2]])
        else:
            pass

index=telabel.index(te)
if len(ref[index][5].split('('))>1:
    rang=int(re.split('[()]',ref[index][5])[1])+int(ref[index][6])
elif len(ref[index][7].split('('))>1:
    rang=int(re.split('[()]',ref[index][7])[1])+int(ref[index][6])

#print(rang)
score_range=[0 for i in range(rang)]#count
control=[0 for i in range(rang)]#control
score_motif=[0 for i in range(rang)]#motif


motif=[]
for i in range(len(motif_a)):
    motif.append([motif_a[i][-5],motif_a[i][-4],motif_a[i][-3],te])

ans=[]
for i in range(len(motif)):
    index=reflabel.index([motif[i][0],motif[i][1],motif[i][2],te])
    print(motif_a[i],te,index,ref[index])
    if float(ref[index][-2])==0.0 and float(ref[index][-3])==0.0:
        ms=int(motif_a[i][1])-int(motif_a[i][-4])
        me=int(motif_a[i][2])-int(motif_a[i][-4])
        if len(ref[index][5].split('('))>1:#reverse
            if  int(ref[index][-5])-1-ms+1>int(ref[index][-5]) or int(ref[index][-5])-1-me+1<int(ref[index][-4]):
                ans.append([str(denovo),int(ref[index][-5])-1-ms+1,int(ref[index][-5])-1-me+1,motif_a[i][-2],motif_a[i][-1],'check'])
            else:
                ans.append([str(denovo),int(ref[index][-5])-1-ms+1,int(ref[index][-5])-1-me+1,motif_a[i][-2],motif_a[i][-1]])
        elif len(ref[index][7].split('('))>1:#forward
            if int(ref[index][-6])-1+ms+1<int(ref[index][-6]) or int(ref[index][-6])+me-1+1>int(ref[index][-5]):
                ans.append([str(denovo),int(ref[index][-6])-1+ms+1,int(ref[index][-6])+me-1+1,motif_a[i][-2],motif_a[i][-1],'check'])
            else:
                ans.append([str(denovo),int(ref[index][-6])-1+ms+1,int(ref[index][-6])+me-1+1,motif_a[i][-2],motif_a[i][-1]])
    else:#insert
        if len(ref[index][5].split('('))>1:
            k=int(ref[index][-1])+4
            length=len([a for a in align[k-4] if a != ''][5])
            bp_del=[]
            while True:
                list_a = [a for a in align[k] if a != '']
                list_g = [a for a in align[k-2] if a != '']
                #print(list_g)
                #print(list_a)
                if list_a[0]=='C':
                    if len(list_g[0].split(':'))>1:
                        now=len(list_g[1])
                        sb=int(list_g[0].split(':')[1].split('-')[0]+'0'*(length-len(list_g[0].split(':')[1].split('-')[0])))+int(list_g[1])
                    else:
                        sb=int(list_g[1])
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
                bp_del_in=bp_del[d]-int(motif_tf[i][1])
                if bp_del_in>=0 and bp_del_in<int(motif_tf[i][2])-int(motif_tf[i][1])+1:
                    b_del.append(bp_del_in)
                elif bp_del_in<0:
                    out_b.append(bp_del_in)
                elif bp_del_in>=int(motif_tf[i][2])-int(motif_tf[i][1])+1:
                    out_o.append(bp_del_in)
            if float(ref[index][-3])!=0.0:#deletion
                k=int(ref[index][-1])+2
                bp_in=[]
                gene_in=[]
                while True:
                    gene_in_b=[]
                    list_g = [a for a in align[k] if a != '']
                    list_a = [a for a in align[k+2] if a != '']
                    #print(list_g)
                    #print(list_a)
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
                        bp_in.append(sb-leave[e]+np.count_nonzero(np.array(enter)<leave[e])-1)
                    if [a for a in align[k+2] if a != ''][-1]==ref[index][7]:
                        break
                    k+=4
                in_b=[]
                in_o=[]
                for m in range(len(gene_in)):
                    if gene_in[m]<int(motif_tf[i][1]):
                        in_b.append(m)
                    elif gene_in[m]>int(motif_tf[i][2])-1:
                        in_o.append(m)
                b_in=len(bp_in)-len(in_b)-len(in_o)
                #print('b_in:',b_in)
                #print('in_b:',in_b)
                #print('in_o:',in_o)
                ms=int(motif_a[i][1])-int(motif_a[i][-4])
                me=int(motif_a[i][2])-int(motif_a[i][-4])
                #print(str(denovo),int(ref[index][-5])-1-ms-len(in_b)+len(out_b)+1,int(ref[index][-5])-1-len(in_b)+len(out_b)+len(b_del)-b_in-me+1,motif_a[i][-2],motif_a[i][-1])
                #print('in_b:',in_b)
                #print('out_b:',out_b)
                #print('b_del:',b_del)
                if int(ref[index][-5])-1-ms-len(in_b)+len(out_b)+1>int(ref[index][-5]) or int(ref[index][-5])-1-len(in_b)+len(out_b)+len(b_del)-b_in-me+1<int(ref[index][-4]):
                    ans.append([str(denovo),int(ref[index][-5])-1-ms-len(in_b)+len(out_b)+1,int(ref[index][-5])-1-len(in_b)+len(out_b)+len(b_del)-b_in-me+1,motif_a[i][-2],motif_a[i][-1],'check'])
                else:
                    ans.append([str(denovo),int(ref[index][-5])-1-ms-len(in_b)+len(out_b)+1,int(ref[index][-5])-1-len(in_b)+len(out_b)+len(b_del)-b_in-me+1,motif_a[i][-2],motif_a[i][-1]])
            else:
                ms=int(motif_a[i][1])-int(motif_a[i][-4])
                me=int(motif_a[i][2])-int(motif_a[i][-4])
                #print(str(denovo),int(ref[index][-5])-ms+len(out_b)-1+1,int(ref[index][-5])-1+len(out_b)+len(b_del)-me+1)
                #print('out_b:',out_b)
                #print('b_del:',b_del)
                if int(ref[index][-5])-ms+len(out_b)-1+1>int(ref[index][-5]) or int(ref[index][-5])-1+len(out_b)+len(b_del)-me+1<int(ref[index][-4]):
                    ans.append([str(denovo),int(ref[index][-5])-ms+len(out_b)-1+1,int(ref[index][-5])-1+len(out_b)+len(b_del)-me+1,motif_a[i][-2],motif_a[i][-1],'check'])
                else:
                    ans.append([str(denovo),int(ref[index][-5])-ms+len(out_b)-1+1,int(ref[index][-5])-1+len(out_b)+len(b_del)-me+1,motif_a[i][-2],motif_a[i][-1]])
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
                if [a for a in align[k] if a != ''][-1]==ref[index][6]:
                    break
                k+=4
            b_del=[]
            out_b=[]
            out_o=[]
            for d in range(len(bp_del)):
                bp_del_in=bp_del[d]-int(motif_tf[i][1])+1
                if bp_del_in>=0 and bp_del_in<int(motif_tf[i][2])-int(motif_tf[i][1])+1:
                    b_del.append(bp_del_in)
                elif bp_del_in<0:
                    out_b.append(bp_del_in)
                elif bp_del_in>=int(motif_tf[i][2])-int(motif_tf[i][1])+1:
                    out_o.append(bp_del_in)
            if float(ref[index][-3])!=0.0:#deletion
                k=int(ref[index][-1])+2
                bp_in=[]
                gene_in=[]
                while True:
                    gene_in_b=[]
                    list_g = [a for a in align[k] if a != '']
                    list_a = [a for a in align[k+2] if a != '']
                    #print(list_g)
                    #print(list_a)
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
                    if gene_in[m]<int(motif_tf[i][1]):
                        in_b.append(m)
                    elif gene_in[m]>int(motif_tf[i][2])-1:
                        in_o.append(m)
                b_in=len(bp_in)-len(in_b)-len(in_o)
                #print('b_in:',b_in)
                #print('in_b:',in_b)
                #print('in_o:',in_o)
                ms=int(motif_a[i][1])-int(motif_a[i][-4])
                me=int(motif_a[i][2])-int(motif_a[i][-4])
                #print(str(denovo),ms+int(ref[index][-6])-1+len(in_b)-len(out_b)+1,me+int(ref[index][-6])-1+len(in_b)-len(out_b)+b_in-len(b_del)+1)
                #print('in_b:',in_b)
                #print('out_b:',out_b)
                #print('b_del:',b_del)
                if ms+int(ref[index][-6])-1+len(in_b)-len(out_b)+1<int(ref[index][-6]) or me+int(ref[index][-6])-1+len(in_b)-len(out_b)+b_in-len(b_del)+1>int(ref[index][-5]):
                    ans.append([str(denovo),ms+int(ref[index][-6])-1+len(in_b)-len(out_b)+1,me+int(ref[index][-6])-1+len(in_b)-len(out_b)+b_in-len(b_del)+1,motif_a[i][-2],motif_a[i][-1],'check'])
                else:
                    ans.append([str(denovo),ms+int(ref[index][-6])-1+len(in_b)-len(out_b)+1,me+int(ref[index][-6])-1+len(in_b)-len(out_b)+b_in-len(b_del)+1,motif_a[i][-2],motif_a[i][-1]])
            else:
                ms=int(motif_a[i][1])-int(motif_a[i][-4])
                me=int(motif_a[i][2])-int(motif_a[i][-4])
                #print(str(denovo),ms+int(ref[index][-6])-len(out_b)-1+1,me+int(ref[index][-6])-1-len(out_b)-len(b_del)-1+1)
                #print('out_b:',out_b)
                #print('b_del:',b_del)
                if ms+int(ref[index][-6])-len(out_b)-1+1<int(ref[index][-6]) or me+int(ref[index][-6])-1-len(out_b)-len(b_del)-1+1>int(ref[index][-5]):
                    ans.append([str(denovo),ms+int(ref[index][-6])-len(out_b)-1+1,me+int(ref[index][-6])-1-len(out_b)-len(b_del)-1+1,motif_a[i][-2],motif_a[i][-1],'check'])
                else:
                    ans.append([str(denovo),ms+int(ref[index][-6])-len(out_b)-1+1,me+int(ref[index][-6])-1-len(out_b)-len(b_del)-1+1,motif_a[i][-2],motif_a[i][-1]])

if denovo!='None':
    with open('out/denovo_'+str(denovo)+'_'+str(te)+'.Freq_refTE_motif_chromVAR.txt','w') as file:
        writer = csv.writer(file,delimiter='\t')
        writer.writerows(ans)
