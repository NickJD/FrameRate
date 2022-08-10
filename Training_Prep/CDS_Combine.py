import gzip
import glob
import collections
import os
import textwrap
import argparse


def fasta_load(sequences,fasta_in):
    first = True
    for line in fasta_in:
        line = line.strip()
        if line.startswith(';'):
            continue
        elif line.startswith('>') and not first:
            sequences.update({sequence_name: seq})
            seq = ''
            sequence_name = line.split()[0]
        elif line.startswith('>'):
            seq = ''
            sequence_name = line.split()[0]
        else:
            seq += str(line)
            first = False
    sequences.update({sequence_name: seq})
    return sequences


sequences = collections.OrderedDict()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='FrameRate: Training Data Preparation.')
    parser.add_argument('-d', action="store", dest='dir_in', required=True,
                        help='Input Directory with file extension separated by comma (./dir_of_cds_files,.fa.gz)')
    parser.add_argument('-lw', action="store", dest='line_wrap', default=True, type=eval, choices=[True, False],
                        help='Default - True: Line wrap FASTA sequence output at 60 chars')
    parser.add_argument('-o', action="store", dest='out_file', required=True,
                        help='Default - Without filetype - default appends \'_StORF-R\' to end of input gff filename (replaces \'.gff\')')
    parser.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')

    options = parser.parse_args()

    if options.gz:
        out_file = gzip.open(options.out_file+'.gz', 'wt', newline='\n', encoding='utf-8')
    else:
        out_file = open(options.out_file, 'w', newline='\n', encoding='utf-8')


    directory_in = options.dir_in.split(',')[0]
    extension_in = options.dir_in.split(',')[1]
    for file_in in glob.glob(directory_in+'/*'+extension_in):
        current_identifier = file_in.split('/')[-1].split(extension_in)[0]
        try:  # Detect whether fasta files are .gz or text and read accordingly
            fasta_in = gzip.open(file_in, 'rt')
        except:
            fasta_in = open(file_in, 'r')
        sequences = collections.OrderedDict()
        sequences = fasta_load(sequences, fasta_in)
        for seq_name, seq in sequences.items():
            seq_name = '>' + current_identifier + ':' + seq_name.replace('>', '')
            out_file.write(seq_name+'\n')
            if options.line_wrap:
                wrapped = textwrap.wrap(seq, width=60)
                for wrap in wrapped:
                    out_file.write(wrap + '\n')
            else:
                out_file.write(seq + '\n')


