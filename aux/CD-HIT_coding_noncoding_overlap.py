import collections



infile = open('/mnt/L_Data/Nextcloud/Dimsum/Results/ensembl_cds_genes_combined_cdhit_eggnogmapped_COGed_6AAs_nostop.filtered','r')

clusters = collections.OrderedDict()

first = True
reps = collections.OrderedDict()
## Load in all data for easier reuse later
for line in infile:
    if line.startswith('>'):
        if first == False:
            cluster_size = len(clusters[cluster_id])

            if len(clusters[cluster_id]) == 1: # Stop at clusters smaller than 10
                break

        cluster_id = line.strip('>')
        cluster_id = cluster_id.strip('\n')
        #cluster_id = cluster_id.split(' ')[1]
        clusters.update({cluster_id: [0,0]})


        first = False
    else:
        clustered = line.split('\t')[1]
        clustered = clustered.split('>')[1]
        clustered = clustered.split('...')[0]
        if '*' in line:
            rep = clustered
            reps.update({rep:[0,0]})
        if first == False:
            if clustered[-1] == '0':
                clusters[cluster_id][1] +=1
            elif clustered[-1] == '1':
                clusters[cluster_id][0] +=1

print("D")

cluster_mixing = list(clusters.values())

coding,non_coding,mixed = 0,0,0
coding_size,non_coding_size,mixed_size = [],[],[]

for cluster in cluster_mixing:
    if cluster[0] >=1 and cluster[1] >=1:
        mixed +=1
        tmp = cluster[0]+cluster[1]
        mixed_size.append(tmp)
    elif cluster[0] >=1 and cluster[1] ==0:
        coding +=1
        coding_size.append(cluster[0])
    if cluster[0] ==0 and cluster[1] >=1:
        non_coding +=1
        non_coding_size.append(cluster[1])

import numpy as np

print(mixed)
print(np.mean(mixed_size))
print(coding)
print(np.mean(coding_size))
print(non_coding)
print(np.mean(non_coding_size))

print(sum(mixed_size))
print(sum(coding_size))
print(sum(non_coding_size))