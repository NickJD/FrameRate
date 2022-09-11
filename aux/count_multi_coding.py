import numpy as np
import collections
import argparse


def count_multi(fasta_in):
    count = collections.defaultdict(int)
    for line in input:
        if line.startswith('>'):
            data = line.split('_')
            if data[0] not in count:
                count[data[0]] = 1
            else:
                count[data[0]] +=1
    counts = list(count.values())
    print(np.mean(counts))
    print(np.median(counts))
    print(counts.count(1))
    print(counts.count(2))
    print(counts.count(3))
    print(counts.count(4))
    print(counts.count(5))
    for key, value in count.items():
        if value == 5:
            print("d")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', action="store", dest='fasta', default="", required=True,
                        help='FASTA to work on')

    options = parser.parse_args()
    input = open(options.fasta,'r')
    count_multi(input)
