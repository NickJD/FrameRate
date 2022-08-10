import gzip
import glob
import collections
import os
import re
import argparse


def revCompIterative(watson):  # Gets Reverse Complement
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                   'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M',
                   'M': 'K', 'V': 'B', 'B': 'V', 'H': 'D', 'D': 'H'}
    watson = watson.upper()
    watsonrev = watson[::-1]
    crick = ""
    for nt in watsonrev:
        try:
            crick += complements[nt]
        except KeyError:
            crick += nt  # Do not modify non-standard DNA
    return crick

gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'*',  ##Stop stars removed
      'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W'}

def check_For_Stops(aa_seq): #Check if frame has stops (*) and if they cut the seq down too far
    stops = []
    stops += [match.start() for match in re.finditer(re.escape('*'), aa_seq)]
    for stop in stops:
        if stop > 5 and stop < 70: # not finished
            return None
    return aa_seq


def translate_frame(sequence):
    translate = ''.join([gencode.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate

def fasta_load(sequences,fasta_in):
    first = True
    for line in fasta_in:
        line = line.strip()
        if line.startswith(';'):
            continue
        elif line.startswith('>') and not first:
            sequences.update({sequence_name: seq})
            if len(sequences) == 100000: # to be removed
                break # and this
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
aa_sequences = collections.OrderedDict()
#out = open('./Extended_CoDing_Sequences_For_Training_Biggest_tmp.faa','w')

count = 0
current_counter = 0
seq_counter = 0

seen_aa = collections.defaultdict()

number_of_CDSs = 0











if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='FrameRate: Training Data Preparation.')
    parser.add_argument('-f', action="store", dest='fasta_in', required=False,
                        help='Input FASTA File')
    parser.add_argument('-d', action="store", dest='dir_in', required=False,
                        help='Input Directory with file extension separated by comma (./dir_of_cds_files,.fa.gz)')
    parser.add_argument('-stop_remove', action="store", dest='stop_remove', default=True, type=eval, choices=[True, False],
                        help='Default - False: Remove stop codons (*) from negative frames ')
    parser.add_argument('-o', action="store", dest='out_file', required=True,
                        help='Default - Without filetype - default appends \'_StORF-R\' to end of input gff filename (replaces \'.gff\')')
    parser.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')

    options = parser.parse_args()

    if options.fasta_in is None and options.dir_in is None:
        parser.error("at least one of -f or -d required")

    if options.gz:
        out_file = gzip.open(options.out_file+'.gz', 'wt', newline='\n', encoding='utf-8')
    else:
        out_file = open(options.out_file, 'w', newline='\n', encoding='utf-8')

    if options.fasta_in:
        try: # Detect whether fasta files are .gz or text and read accordingly
            fasta_in = gzip.open(options.fasta_in,'rt')
        except:
            fasta_in = open(options.fasta_in,'r')
        sequences = collections.OrderedDict()
        sequences = fasta_load(sequences, fasta_in,'combined')
    elif options.dir_in:
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
            number_of_CDSs += len(sequences)
            for seq_name, seq in sequences.items():
                seq_name = seq_name.replace('>', '')
                rev_seq = revCompIterative(seq)
                ### Frame 1
                aa_seq = translate_frame(seq)
                aa_seq = aa_seq[1:]
                if aa_seq not in seen_aa:
                    seen_aa[aa_seq] = None
                    out_file.write('>' + current_identifier + '_' + seq_name + '_' + str(seq_counter) + ',' + str(count) + ',1\n' + aa_seq + '\n')
                    seq_counter += 1
                ### Frame 2
                tmp_seq = seq[1:]
                aa_seq = translate_frame(tmp_seq)
                aa_seq = aa_seq[1:]
                if aa_seq not in seen_aa:
                    seen_aa[aa_seq] = None
                    out_file.write('>' + current_identifier + '_' + seq_name + '_' + str(seq_counter) + ',' + str(count) + ',0\n' + aa_seq + '\n')
                    seq_counter += 1
                ### Frame 3
                tmp_seq = tmp_seq[1:]
                aa_seq = translate_frame(tmp_seq)
                aa_seq = aa_seq[1:]
                if aa_seq not in seen_aa:
                    seen_aa[aa_seq] = None
                    out_file.write('>' + current_identifier + '_' + seq_name + '_' + str(seq_counter) + ',' + str(count) + ',0\n' + aa_seq + '\n')
                    seq_counter += 1
                ### Frame 4
                aa_seq = translate_frame(rev_seq)
                aa_seq = aa_seq[1:]
                if aa_seq not in seen_aa:
                    seen_aa[aa_seq] = None
                    out_file.write('>' + current_identifier + '_' + seq_name + '_' + str(seq_counter) + ',' + str(count) + ',0\n' + aa_seq + '\n')
                    seq_counter += 1
                ### Frame 5
                tmp_rev_seq = rev_seq[1:]
                aa_seq = translate_frame(tmp_rev_seq)
                aa_seq = aa_seq[1:]
                if aa_seq not in seen_aa:
                    seen_aa[aa_seq] = None
                    out_file.write('>' + current_identifier + '_' + seq_name + '_' + str(seq_counter) + ',' + str(count) + ',0\n' + aa_seq + '\n')
                    seq_counter += 1
                ### Frame 6
                tmp_rev_seq = tmp_rev_seq[1:]
                aa_seq = translate_frame(tmp_rev_seq)
                aa_seq = aa_seq[1:]
                if aa_seq not in seen_aa:
                    seen_aa[aa_seq] = None
                    out_file.write('>' + current_identifier + '_' + seq_name + '_' + str(seq_counter) + ',' + str(count) + ',0\n' + aa_seq + '\n')
                    seq_counter += 1

                count += 1
            current_counter += 1
            break
