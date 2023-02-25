import csv
import numpy as np


with open('out/denovo_1_motif_pos_TEsubfamily.txt') as f: 
    reader = csv.reader(f, delimiter='\t')
    TE_peak = [row for row in reader]#TSS_peak.generate.txt

te='MER130'

ans=[]
for i in range(len(TE_peak)):
    if TE_peak[i][-2]==te:
        ans.append([TE_peak[i][0],TE_peak[i][1],TE_peak[i][2]])

with open('out/denovo_1_motif_pos_'+str(te)+'.bed','w') as file:
    writer = csv.writer(file,delimiter='\t')
    writer.writerows(ans)





with open('out/denovo_3_motif_pos_TEsubfamily.txt') as f: 
    reader = csv.reader(f, delimiter='\t')
    TE_peak = [row for row in reader]#TSS_peak.generate.txt

te='MamRep434'

ans=[]
for i in range(len(TE_peak)):
    if TE_peak[i][-2]==te:
        ans.append([TE_peak[i][0],TE_peak[i][1],TE_peak[i][2]])

with open('out/denovo_3_motif_pos_'+str(te)+'.bed','w') as file:
    writer = csv.writer(file,delimiter='\t')
    writer.writerows(ans)