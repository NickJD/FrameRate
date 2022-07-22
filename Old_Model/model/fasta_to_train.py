

input = open('./Extended_CoDing_Sequences_For_Training_Biggest.faa_CD_60.fa','r')
out = open('./Extended_CoDing_Sequences_For_Training_Biggest.faa_CD_60.csv','w')

first = True
for line in input:
    line = line.strip()

    if line.startswith('>') and not first:
        out.write(sequence_name+','+seq+'\n')
        seq = ''
        line = line.replace('>', '')
        sequence_name = line
    elif line.startswith('>'):
        seq = ''
        line = line.replace('>', '')
        sequence_name = line
    else:
        seq += str(line)
        first = False

