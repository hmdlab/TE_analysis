#!/usr/bin/env python
# coding: utf-8

import csv
import numpy as np



import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-tf','--tf',type=str,default='None',help='TFname')
parser.add_argument('-te','--te',type=str,default='None',help='TEname')
parser.add_argument('-tf_m','--tf_m',type=str,default='None',help='Selecet denovo motif number like denovo_1')
args=parser.parse_args()
tf=args.tf 
te=args.te 
tf_m=args.tf_m



#peak
f = open('out/'+str(tf)+'_'+str(te)+'.Freq_refTE_count.txt')
data3 = f.read() 
f.close()
score_rang= data3.split('\n')

score_rang=np.array(score_rang)
score_rang=np.delete(score_rang,-1)
score_rang=score_rang.astype('float32')

f = open('out/'+str(tf)+'_'+str(te)+'.Freq_refTE_control.txt')
data3 = f.read() 
f.close()
control= data3.split('\n')

control=np.delete(control,-1)
control=control.astype('float32')


with open('out/'+str(tf_m)+'_'+str(te)+'.Freq_refTE_motif_chromVAR.txt') as f:
    reader = csv.reader(f, delimiter='\t')
    motif = [row for row in reader]




for i in range(len(motif)):
    motif[i][1]=int(motif[i][1])
    motif[i][2]=int(motif[i][2])


f = open('out/'+str(tf)+'_'+str(te)+'.phyloPsocre_refTE_motif_mapping.txt')
data3 = f.read()
f.close()
motif= data3.split('\n')

for i in range(len(score_rang)):
    if control[i]!=0.0:
        control[i]=score_rang[i]/control[i]


#Freq
import numpy as np
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 22


motif_c=[['',216,223]]#MER130-Neurod2
#motif_c=[['',222,229]]#Mamrep434-Lhx2

left_1=[]
for i in range(len(motif_c)):
    for j in range(max(int(motif_c[i][1]),int(motif_c[i][2]))-min(int(motif_c[i][1]),int(motif_c[i][2]))+1):
        left_1.append(min(int(motif_c[i][1]),int(motif_c[i][2]))+j)

left_1=list(set(left_1))
left_1=np.sort(left_1)
left_1=[s for s in left_1 if 0<s and s<=len(control)]
height_m=[]
left_m=[]
for i in range(len(left_1)-1):
    height_m.append(height[int(left_1[i])-1])
    left_m.append(int(left_1[i]))
    if left_1[i+1]==None:
        break
    if int(left_1[i+1])-int(left_1[i])!=1:
        left_m.append(None)
        height_m.append(None)


for i in range(len(height_m)):
    height_m[i]=0



l=[]
for i in range(len(score_rang)):
    l.append(i+1)
left = np.array(l)
height_1 = np.array(score_rang)
plt.figure(figsize=(15,5),dpi=400)
plt.xlabel(str(te)+'   [bp]')
plt.ylabel('ChIP-seq read count')#coverage
p1 = plt.plot(left, height_1)
#plt.savefig("out/"+str(te)+"_ChIP_ref.png", bbox_inches="tight")
plt.show()





# normalize

height_m=[]
left_m=[]
for i in range(len(left_1)-1):
    height_m.append(height_n[int(left_1[i])-1])
    left_m.append(int(left_1[i]))
    if left_1[i+1]==None:
        break
    if int(left_1[i+1])-int(left_1[i])!=1:
        left_m.append(None)
        height_m.append(None)



for i in range(len(height_m)):
    height_m[i]=0


l=[]
for i in range(len(score_rang)):
    l.append(i+1)
left = np.array(l)
height_n = np.array(control*100)
plt.figure(figsize=(15,5),dpi=400)
plt.xlabel(str(te)+'   [bp]')
plt.ylabel('ChIP-seq peak count /\nTE fragment count\n[%]')
p1 = plt.plot(left, height_n)
#plt.savefig("out/"+str(te)+"_ChIP_ref_norm.png", bbox_inches="tight")
plt.show()


   
fig = plt.figure(figsize=(16,14),dpi=400)

ax1 = fig.add_subplot(5,1,(1,2))
l=[]
for i in range(len(score_rang)):
    l.append(i+1)
left = np.array(l)
height_1 = np.array(score_rang)
ax1.set_ylabel('ChIP-seq peak\ncount')
ax1.plot(left, height_1,label=str(tf))
#ax1.bar(left_m, [max(height_1)+1]*len(height_m),label='$\it{de}$ $\it{novo}$ ' +str(tf_m.split('_')[1]),width=1,color='tab:orange')
#ax1.bar(left_m, [-0.5]*len(height_m),width=1,color='tab:orange')
ax1.set_ylim(-0.2, max(height_1)+0.2)
#ax1.legend(loc=2)
ax1.set_title(str(tf)+' ('+str(tf_m)+')')
ax1.set_title(str(tf)+' ($\it{de}$ $\it{novo}$ ' +str(tf_m.split('_')[1])+')')


ax2 = fig.add_subplot(5,1,(3,4))
height_n = np.array(control*100)
ax2.set_ylabel('ChIP-seq peak count /\nTE fragment count\n[%]')
ax2.plot(left, height_n,label=str(tf))
#ax2.bar(left_m, [max(height_n)+1]*len(height_m),label='$\it{de}$ $\it{novo}$ ' +str(tf_m.split('_')[1]),width=1,color='tab:orange')
#ax2.bar(left_m, [-0.5]*len(height_m),width=1,color='tab:orange')
ax2.set_ylim(-0.2, max(height_n)+0.5)
#ax2.legend(loc=2)

ax3 = fig.add_subplot(5,1,5)
left_1=[]
for i in range(len(motif)):
    for j in range(max(int(motif[i][1]),int(motif[i][2]))-min(int(motif[i][1]),int(motif[i][2]))+1):
        left_1.append(min(int(motif[i][1]),int(motif[i][2]))+j)

height_2=[0]*len(left)
for i in range(len(left)):
    height_2[i]=left_1.count(i)
height_2m=[0]*len(left)
for i in range(motif_c[0][1],motif_c[0][2]+1):
    height_2m[i]=left_1.count(i)

ax3.bar(left, height_2,width=1)
ax3.set_xlabel(str(te)+'   [bp]')
ax3.set_ylabel('accessible\nmotif count')
#plt.savefig("out/"+str(tf)+"_"+str(tf_m)+"_"+str(te)+"_motif_ref.png", bbox_inches="tight")
plt.show()

