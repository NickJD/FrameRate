

input = open('./Extended_CoDing_Sequences_For_Training_CD_80.fa','r')
out = open('./Extended_CoDing_Sequences_For_Training.csv','w')

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

