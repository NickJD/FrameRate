import argparse
import collections
import keras
import tensorflow.compat.v1 as tf
import re
from utils import *
from itertools import islice

###################
gencode = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}

aaletters = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
aaindex = {'K': 0, 'H': 1, 'T': 2, 'V': 3, 'Q': 4, 'L': 5, 'D': 6, 'R': 7, 'A': 8, 'G': 9, 'Y': 10, 'F': 11, 'S': 12,
           'N': 13, 'P': 14, 'C': 15, 'M': 16, 'I': 17, 'W': 18, 'E': 19}

complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


# , 'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M','M': 'K', 'V': 'B', 'B': 'V', 'H': 'D', 'D': 'H'}
############################

def convert_to_frames(seq_id, seq, reads):
    global count_none
    none_flag = True
    for i in range(3):
        aa_seq = translate_frame(seq[i:])
        aa_seq_ln = longest_nostop(aa_seq)

        if aa_seq_ln is not None:
            reads[seq_id].append([seq_id + '_frame_' + str(i), aa_seq_ln])
            none_flag = False
    for i in [2, 1, 0]:  # reverse compliments
        rev_seq = revCompIterative(seq[:len(seq) - i])
        aa_seq2 = translate_frame(rev_seq)
        aa_seq2_ln = longest_nostop(aa_seq2)

        if aa_seq2_ln is not None:
            reads[seq_id].append([seq_id + '_frame_' + str(6 - i), aa_seq2_ln])
            none_flag = False

    if none_flag == True:
        count_none += 1


def translate_frame(sequence):
    translate = ''.join([gencode[sequence[3 * i: 3 * i + 3]] for i in range(len(sequence) // 3)])
    return translate


def longest_nostop(translate):
    nostop_region_length, nostop_region = max(
        (len(ss), ss) for ss in translate.split('*'))  # Get longest segment of AA sequence without stop codons
    if nostop_region_length >= options.min_frame:
        nostop_region = nostop_region[-75:]  # Get the last 75 AAs - Most likely not needed
        return nostop_region
    else:
        return None


def revCompIterative(watson):
    watson = watson.upper()
    watsonrev = watson[::-1]
    crick = ""
    for nt in watsonrev:
        crick += complements[nt]
    return crick


def predictor(reads, model, options, coding, noncoding):
    batch_size = options.batch_size

    ids = []
    x_pred = []
    for read, frames in reads.items():
        for frame in frames:
            ids.append(frame[0])
            x_pred.append(to_onehot(frame[1], aaindex))

    print("The number of reads to be classified is", len(ids))

    x_pred = np.array(x_pred)
    print("The shape of input data is", x_pred.shape)
    y_pred = model.predict(x_pred, verbose=1, batch_size=batch_size)

    print("Prediction done, start writing output files")

    for i, key in enumerate(ids):
        current_read = key.split('_')[0]
        frames = reads[current_read]
        for j, frame in enumerate(frames):
            if key in frame:
                reads[current_read][j].append(y_pred[i][0])

    for read, frames in reads.items():
        frames_sorted = sorted(frames, key=lambda x: x[2], reverse=True)  # Frame with highest score
        for frame in frames_sorted:
            if frame[2] < 0.5:
                noncoding.write(frame[0] + '_Score:' + str(format(frame[2], '.2f')) + '\n' + frame[1] + '\n')
            else:
                coding.write(frame[0] + '_Score:' + str(format(frame[2], '.2f')) + '\n' + frame[1] + '\n')

    print("Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', action="store", dest='fasta', default="", required=True,
                        help='FASTA to work on')
    parser.add_argument('-c', action="store", dest='chunk', default=100000, type=int, required=False,
                        help='Number of Reads to for each Chunk')
    parser.add_argument('-m', '--model_file', action='store', dest='model_file', required=False,
                        help='Pretrained model to use')
    parser.add_argument('-min_frame', action="store", dest='min_frame', default=50, type=int, required=False,
                        help='Default - 50: Minimum frame size in AA')
    parser.add_argument('-o', '--output_prefix', action='store', dest='out_prefix', required=True,
                        help='Output file prefix')
    parser.add_argument('-b', '--batch_size', action='store', dest='batch_size', type=int, default=5000, required=False,
                        help='Batch size for the model')
    parser.add_argument('-s', '--subsample', action='store', dest='subsample', type=float, default=1.0, required=False,
                        help='Subsample factor for the input file')

    options = parser.parse_args()

    print("Loading saved model", options.model_file)
    model = keras.models.load_model(options.model_file)

    print("Starting reading input ", options.fasta)
    coding = open(options.out_prefix + '_coding.fa', 'w')
    noncoding = open(options.out_prefix + '_noncoding.fa', 'w')

    piece = 0
    piece_size = options.chunk * 2
    with open(options.fasta, 'r') as f:
        for lines in iter(lambda: list(islice(f, piece_size)), []):
            piece += 1

            if np.random.uniform() <= options.subsample:
                print("Processing file piece", piece, "with", options.chunk, "reads")

                reads = collections.defaultdict(list)
                count, count_none = 0, 0
                for i in range(len(lines)):
                    line = lines[i].strip()
                    if line.startswith('>'):
                        sequence_name = line
                        seq = lines[i + 1].strip()
                        try:
                            convert_to_frames(sequence_name, seq, reads)
                        except KeyError:
                            count += 1

                print("The number of valid sequences after filtering in the input file is ", len(reads))
                print("The number of reads that have non-canonical bases is ", count)
                print("The number of genes that did not pass the minimum frame length filtering is ", count_none)

                n_filtered_frames = 0
                for frames in reads.values():
                    n_filtered_frames += len(frames)
                print("The number of frames is ", n_filtered_frames)
                #print("The frame to sequence ratio is ", n_filtered_frames / len(reads))

                predictor(reads, model, options, coding, noncoding)

            else:
                print("Skipping file piece", piece)
    coding.close()
    noncoding.close()