import collections
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd

CDS_Coding = open('/home/nick/Git/CoDing_Sequence_Frame_Prediction_From_DNA_Reads/classifier/Running'
             '/Coding.fa','r')
CDS_Non_Coding = open('/home/nick/Git/CoDing_Sequence_Frame_Prediction_From_DNA_Reads/classifier/Running'
             '/Non_Coding.fa','r')

Unassembled_Coding = open('/home/nick/Git/CoDing_Sequence_Frame_Prediction_From_DNA_Reads/classifier/Finished'
             '/megahit_assembly_contigs_Min_1000_aligned_unaligned_reads_Coding.fa','r')
Unassembled_Non_Coding = open('/home/nick/Git/CoDing_Sequence_Frame_Prediction_From_DNA_Reads/classifier/Finished'
             '/megahit_assembly_contigs_Min_1000_aligned_unaligned_reads_Non_Coding.fa','r')

CDS_Coding_scores = []
for line in CDS_Coding:
    if line.startswith('>'):
        score = line.split(':')[-1].strip()
        CDS_Coding_scores.append(float(score))
CDS_Coding_scores = [ round(elem, 2) for elem in CDS_Coding_scores]
print(np.mean(CDS_Coding_scores))
print(np.median(CDS_Coding_scores))
######################
# CDS_Non_Coding_scores = []
# for line in CDS_Non_Coding:
#     if line.startswith('>'):
#         score = line.split(':')[-1].strip()
#         CDS_Non_Coding_scores.append(float(score))
# print(np.mean(CDS_Non_Coding_scores))
# print(np.median(CDS_Non_Coding_scores))
# ######################
# Unassembled_Coding_Scores = []
# for line in Unassembled_Coding:
#     if line.startswith('>'):
#         score = line.split(':')[-1].strip()
#         Unassembled_Coding_Scores.append(float(score))
# print(np.mean(Unassembled_Coding_Scores))
# print(np.median(Unassembled_Coding_Scores))
# #######################
# Unassembled_Non_Coding_Scores = []
# for line in Unassembled_Non_Coding:
#     if line.startswith('>'):
#         score = line.split(':')[-1].strip()
#         Unassembled_Non_Coding_Scores.append(float(score))
# print(np.mean(Unassembled_Non_Coding_Scores))
# print(np.median(Unassembled_Non_Coding_Scores))

CDS_Coding_scores_DF = pd.DataFrame({'CDS_Coding_scores':pd.Series(CDS_Coding_scores)})
# CDS_Non_Coding_scores_DF = pd.DataFrame({'CDS_Non_Coding_scores':pd.Series(CDS_Non_Coding_scores)})
# Unassembled_Coding_DF = pd.DataFrame({'Unassembled_Coding':pd.Series(Unassembled_Coding_Scores)})
# Unassembled_Non_Coding_DF = pd.DataFrame({'Unassembled_Non_Coding':pd.Series(Unassembled_Non_Coding_Scores)})

CDS_Coding_scores_DF = pd.melt(CDS_Coding_scores_DF)
CDS_Coding_scores_DF.insert(0,'Group','CDS Aligned Coding')

#data = [CDS_Coding_scores_DF,CDS_Non_Coding_scores_DF,Unassembled_Coding_DF,Unassembled_Non_Coding_DF]

#All_Data = pd.concat(data)

#ax = sns.boxplot(data=data)
plt.show()