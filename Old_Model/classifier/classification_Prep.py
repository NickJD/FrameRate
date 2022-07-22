import re
import collections

def revCompIterative(watson): #Gets Reverse Complement
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                   'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M',
                   'M': 'K', 'V': 'B', 'B': 'V', 'H': 'D', 'D': 'H'}
    watson = watson.upper()
    watsonrev = watson[::-1]
    crick = ""
    for nt in watsonrev: # make dict to catch bad nts - if more than 1 output them to std error
        try:
            crick += complements[nt]
        except KeyError:
            crick += nt
    return crick


def find_stops(seq):
    stops = []
    for stop_codon in "TAG,TGA,TAA".split(','):  # Find all Stops in seq
        stops += [match.start() for match in re.finditer(re.escape(stop_codon), seq)]
    return stops



input = open('Escherichia_coli_Test_Data.fa','r')

sequences = collections.OrderedDict()

first = True
for line in input:
    line = line.strip()
    if line.startswith('>') and not first:
        sequences.update({sequence_name: seq})
        seq = ''
        sequence_name = line
    elif line.startswith('>'):
        seq = ''
        sequence_name = line
    else:
        seq += str(line)
        first = False
sequences.update({sequence_name: seq})


for sequence_name, seq in sequences.items():
    fw_stops = find_stops(seq)
    seq_rev = revCompIterative(seq)
    bw_stops = find_stops(seq_rev)
    print("D")
    #If found stops make seq less than 70aa - is 70, make up t o 75 with first 5 of seq.
    #Trim seqs here at 50 aa











