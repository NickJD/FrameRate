import collections
import sys
import operator
import tqdm


seqs = []

seqs_front = collections.defaultdict(int)
seqs_rear = collections.defaultdict(int)
seq_IDs = collections.defaultdict(list)

fasta_in = open('/home/nick/Documents/tmp/out_retrained_model_coding.fa','r')

count = 1
for line in fasta_in:  # .split('\n'):
    line = line.strip()
    if line.startswith('>'):# and first == False:  # Check if first seq in file
        ID = line.strip()
        count +=1
    else:
        seq = line.strip()
        seqs.append(seq)
        #seqs_front[seq[0:25]] = count
        #seqs_rear[seq[-25:]] = count
        if seq not in seq_IDs:
            seq_IDs[seq] = [ID,seq,0]
        else:
            seq_IDs[seq][2] +=1

#######################################################




######################

assemblies = collections.defaultdict(list)

already_assembled = collections.defaultdict()

for seq, data  in tqdm.tqdm(seq_IDs.items()):
    if data[0] not in already_assembled:
        front_matches = [seq_IDs[key] for key in seq_IDs.keys() if seq[-25:] in key[0:25]]
        if len(front_matches) >=2: # first check if one is a substring of the other.
            frame_sequences = [el[1] for el in front_matches]
            if len({len(i) for i in frame_sequences}) == 1: # take the first - not ideal
                assembled = seq + front_matches[0][1][25:]
                already_assembled[data[0]] = None
                already_assembled[front_matches[0][0]] = None
            else:
                largest_frame = max(front_matches, key=operator.itemgetter(1))[0]
                largest_match = next((l for l in front_matches if largest_frame in l), None)
                assembled = seq + largest_match[1][25:]
                already_assembled[data[0]] = None
                already_assembled[largest_frame] = None
            assemblies[data[0]] = [assembled]
        if front_matches:
            assembled = seq + front_matches[0][1][25:]
            already_assembled[data[0]] = None
            already_assembled[front_matches[0][0]] = None
            assemblies[data[0]] = [assembled]



#########################################
print(len(assemblies))
for key, value in assemblies.items():
    print(key+'\n'+value[0])



