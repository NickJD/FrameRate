import collections


def all_same(items):
    return all(x == items[0] for x in items)


PEP_In = open('Extended_CoDing_Sequences_For_Training_Biggest.faa_CD.clstr','r')

clustered_types = collections.defaultdict(list)
first = True


for line in PEP_In:
    if line.startswith('>'):
        if first == False:

            if len(clustered_types[cluster_id]) == 1: # Stop at clusters smaller than 10
                break

        cluster_id = line.strip('>')
        cluster_id = cluster_id.strip('\n')
        cluster_id = cluster_id.split(' ')[1]
        first = False
    else:
        clustered = line.split('\t')[1]
        clustered = clustered.split('...')[0]
        clustered = clustered.split(',')[-1]

        if first == False:
            clustered_types[cluster_id].append(clustered)

print("Done")


for  key, value in clustered_types.items():
    if all_same(value):
        continue
    else:
        print(key)


