#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import csv
import pandas as pd
import numpy as np

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-denovo','--denovo',type=str,default='None',help='number of de novo motif')
args=parser.parse_args()
denovo=args.denovo 

with open('denovo_'+str(denovo)+'_motif_pos.txt') as f:
    reader = csv.reader(f, delimiter=' ')
    var_all = [row for row in reader]


ans=[]
for i in range(1,len(var_all)):
    ans.append([var_all[i][1],var_all[i][2],var_all[i][3],var_all[i][-2],var_all[i][-1]])

with open('denovo_'+str(denovo)+'_motif_pos.txt','w') as file:
    writer = csv.writer(file,delimiter='\t')
    writer.writerows(ans)

