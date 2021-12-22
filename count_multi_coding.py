import numpy as np
import collections

input = open('/home/nick/Git/CoDing_Sequence_Frame_Prediction_From_DNA_Reads/classifier/Finished/trimmed_paired_SRR873595_combined_subsampled_10_Coding.fa','r')

count = collections.defaultdict(int)

last = ""

for line in input:
#    if ">SRR873595:::126731668:0:0:0:;0.999093" in line:
    if line.startswith('>'):
        data = line.split('_')
        if data[0] not in count:
            count[data[0]] = 1
        else:
            count[data[0]] +=1
        last = data[0]

counts = list(count.values())
print(np.mean(counts))
print(np.median(counts))


