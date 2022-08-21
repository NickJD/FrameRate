import argparse
import collections
import keras
import tensorflow.compat.v1 as tf
from predictor import predictor
import re
from utils import *

def remove_char(str, n):
    first_part = str[:n]
    last_part = str[n + 1:]
    return first_part + last_part

def check_For_Stops(aa_seq): #Check if frame has stops (*) and if they cut the seq down too far
    unstopped_region_length, unstopped_region = max((len(ss), ss) for ss in aa_seq.split('*')) # Get longest segment of AA sequence without stop codons
    unstopped_region = unstopped_region[-75:] # Get the last 75 AAs - Most likely not needed
    if unstopped_region_length >= options.min_frame:
        return unstopped_region
    else:
        return None


def convert_To_Frames(seq_id,seq,Reads):
    #seq_id = seq_id.replace('>', '')
    #rev_seq = revCompIterative(seq)

    aa_seq = translate_frame(seq)
    aa_seq = check_For_Stops(aa_seq)
    if aa_seq != None:
        Reads[seq_id].append([seq_id + '_Frame:1', aa_seq])

    ######################
    aa_seq = translate_frame(seq[1:])
    aa_seq = check_For_Stops(aa_seq)
    if aa_seq != None:
        Reads[seq_id].append([seq_id + '_Frame:2', aa_seq])

    #######################
    aa_seq = translate_frame(seq[2:])
    aa_seq = check_For_Stops(aa_seq)
    if aa_seq != None:
        Reads[seq_id].append([seq_id + '_Frame:3', aa_seq])

    ##################################################
    aa_seq = translate_frame(revCompIterative(seq))
    aa_seq = check_For_Stops(aa_seq)
    if aa_seq != None:
        Reads[seq_id].append([seq_id + '_Frame:4', aa_seq])

    #################################
    aa_seq = translate_frame(revCompIterative(seq[:len(seq)-1]))
    aa_seq = check_For_Stops(aa_seq)
    if aa_seq != None:
        Reads[seq_id].append([seq_id + '_Frame:5', aa_seq])

    #################################
    aa_seq = translate_frame(revCompIterative(seq[:len(seq)-2]))
    aa_seq = check_For_Stops(aa_seq)
    if aa_seq != None:
        Reads[seq_id].append([seq_id + '_Frame:6', aa_seq])

    return Reads

###################
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
      'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
      'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
############################

def translate_frame(sequence):
    translate = ''.join([gencode.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate

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


def DNA_To_Frames(fasta_in):#,chunk_line):

    count = 0
    first = True
    Reads = collections.defaultdict(list)
    #line = linecache.getline(fasta_in, chunk_line)
    for line in fasta_in:#.split('\n'):
        line = line.strip()
        if line.startswith('>') and first == False:  # Check if first seq in file
            Reads = convert_To_Frames(seq_id,seq,Reads)
            count +=1
            seq = ''
            seq_id = line
            #if len(Reads) == int(chunk_size):
            #     print("CHUNK: " + str(chunk))
            #     #predictor(Reads, model)
            #     #del Reads
            #     #Reads = collections.defaultdict(list)
            #     chunk +=1
                #Reads = convert_To_Frames(seq_id, seq, Reads)
            #    break
        elif line.startswith('>'):
            seq = ''
            seq_id = line
        else:
            seq += str(line)
            first = False

    #predictor(Reads, model)
    Reads = convert_To_Frames(seq_id, seq, Reads)


    Filtered_Reads = len(Reads)
    Filtered_Frames = 0
    for listElem in list(Reads.values()):
        Filtered_Frames += len(listElem)
    per_read = (Filtered_Frames/Filtered_Reads)
    print(count)
    print(Filtered_Reads)
    print(Filtered_Frames)
    print(per_read)

    return Reads


def max_value(inputlist):
    return max([sublist[-1] for sublist in inputlist])




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', action="store", dest='fasta', default="", required=True,
                        help='FASTA to work on')
    parser.add_argument('-c', action="store", dest='chunk', default=100000, required=False,
                        help='Number of Reads to for each Chunk')
    parser.add_argument('-m', '--model_file', action='store', dest='model_file', required=False,
                        help='Pretrained model to use')
    # parser.add_argument('-in_stops', action="store", dest='in_stops', default=False, type=eval, choices=[True, False],
    #                     help='Default - False: Remove frames with Stops codons between positions 10-65 ')
    # parser.add_argument('-stops', action="store", dest='stops', default=False, type=eval, choices=[True, False],
    #                     help='Default - True: Extract longest region without Stops codons')
    parser.add_argument('-min_frame', action="store", dest='min_frame', default=50, type=int,
                        help='Default - 50: Minimum frame size in AA')
    parser.add_argument('-o', '--output_prefix', action='store', dest='out_prefix',
                        help='Output file prefix')
    parser.add_argument('-GPU', action="store", dest='gpu', default=False, type=eval, choices=[True, False],
                        help='Default - False: Use GPU for computation')
    options = parser.parse_args()

    if options.gpu == True:
        tf.config.set_visible_devices([], 'GPU') # not implemented yet and will fail if user has not setup tensorflow GPU


    ##### Convert DNA reads into AA frames
    #fasta_in = open(options.fasta, 'r')
    # with open(options.fasta, "r") as f:
    #     while chunk := f.read(int(options.chunk)):  # you can use any chunk size you want
    #         #
    #         Reads = DNA_To_Frames(model,chunk)
    #         print("W")

    with open(options.fasta, "r") as fasta_in:
        current_chunk_num = 0
        lines = []
        first = True
        for line in fasta_in:
            line = line.strip()
            if line.startswith(';'):
                continue
            elif line.startswith('>') and not first:
                lines.append(sequence_name)
                lines.append(seq)
                seq = ''
                sequence_name = line
            elif line.startswith('>'):
                seq = ''
                sequence_name = line
            else:
                seq += str(line)
                first = False
            current_chunk_num +=1

        lines.append(line)
        Reads = DNA_To_Frames(lines)
        print("Reads")
        model = keras.models.load_model(options.model_file)
        predictor(Reads, model, options)

    # try: # Detect whether fasta/gff files are .gz or text and read accordingly
    #     fasta_in = gzip.open(options.fasta,'rt')
    #     DNA_To_Frames(options.chunk,model,fasta_in)
    # except:
    #     fasta_in = open(options.fasta,'r')
    #     DNA_To_Frames(options.chunk,model,fasta_in)



    print("Done Done")





