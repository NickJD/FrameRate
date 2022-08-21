import argparse
import collections
import gzip

def fasta_load(sequences,fasta_in):
    first = True
    for line in fasta_in:
        line = line.strip()
        if line.startswith(';'):
            continue
        elif line.startswith('>') and not first:
            sequences.update({sequence_name: [seq,['']]})
            if len(sequences) == 100000: # to be removed
                break # and this
            seq = ''
            line = line.replace('>', '')
            sequence_name = line.split()[0]
        elif line.startswith('>'):
            seq = ''
            line = line.replace('>', '')
            sequence_name = line.split()[0]
        else:
            seq += str(line)
            first = False
    sequences.update({sequence_name: [seq,['']]})
    return sequences

def COG_load(cog_in,sequences):
    for line in cog_in:
        if not line.startswith('#'):
            line_data = line.split('\t')
            if line_data[0] in sequences:
                COG_Family = line_data[4].split('@')[0]
                sequences[line_data[0]][1] = COG_Family


    return sequences

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', action='store', dest='fasta_in', required=True,
                        help='Original Ensembl CDS FASTA')
    parser.add_argument('-c', action='store', dest='cog_in', required=False,
                        help='Coding EggNog Annotations')
    parser.add_argument('-o', action="store", dest='out_file', required=False,
                        help='Prefix output filename to be used to output FASTAs with and without COGs')
    parser.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')
    options = parser.parse_args()

    try:  # Detect whether fasta files are .gz or text and read accordingly
        fasta_in = gzip.open(options.fasta_in, 'rt')
    except:
        fasta_in = open(options.fasta_in, 'r')
    sequences = collections.OrderedDict()
    sequences = fasta_load(sequences, fasta_in)


    cog_in = open(options.cog_in, 'r')

    Sequences = COG_load(cog_in, sequences)


    if options.gz:
        cogged_out_file = gzip.open(options.out_file+'_COGed.fasta.gz', 'wt', newline='\n', encoding='utf-8')
        not_cogged_out_file = gzip.open(options.out_file + '_Not_COGed.fasta.gz', 'wt', newline='\n', encoding='utf-8')
    else:
        cogged_out_file = open(options.out_file + '_COGed.fasta', 'w', newline='\n', encoding='utf-8')
        not_cogged_out_file = open(options.out_file+'_Not_COGed.fasta', 'w', newline='\n', encoding='utf-8')

    print("w")
    for ID, data in sequences.items():
        if type(data[1]) is str:
            cogged_out_file.write('>' + data[1] + '_' + ID + '\n' + data[0] + '\n')
        else:
            not_cogged_out_file.write('>' + 'NoCOG_' + ID + '\n' + data[0] + '\n')



