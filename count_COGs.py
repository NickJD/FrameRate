import argparse
import collections

from scipy.stats import chisquare



def anno_load(anno_in,COGs):
    count = 0
    for line in anno_in:
        if line.startswith('#'):
            continue
        else:
            line = line.split('\t')
            COG = line[6]
            COGs[COG]+=1
            count +=1

    return COGs,count








if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', action='store', dest='genes', required=True,
                        help='Gene EggNog Annotations')
    parser.add_argument('-c', action='store', dest='coding', required=True,
                        help='Coding EggNog Annotations')
    parser.add_argument('-nc', action='store', dest='non_coding', required=True,
                        help='Non_Coding EggNog Annotations')
    options = parser.parse_args()

    Gene_COGs = collections.defaultdict(int)
    Gene_COGs_P = collections.defaultdict()
    Coding_COGs = collections.defaultdict(int)
    Coding_COGs_P = collections.defaultdict()
    Non_Coding_COGs = collections.defaultdict(int)
    Non_Coding_COGs_P = collections.defaultdict()

    gene_anno_in = open(options.genes, 'r')
    coding_anno_in = open(options.coding, 'r')
    non_coding_anno_in = open(options.non_coding, 'r')
    Gene_COGs, num_Gene_COGs = anno_load(gene_anno_in, Gene_COGs)
    Coding_COGs,num_Coding_COGs = anno_load(coding_anno_in,Coding_COGs)
    Non_Coding_COGs,num_Non_Coding_COGs = anno_load(non_coding_anno_in, Non_Coding_COGs)

    for COG, num in Gene_COGs.items():
        Gene_COGs_P[COG] = float("{:.2f}".format((num/num_Gene_COGs) * 100))
    for COG, num in Coding_COGs.items():
        Coding_COGs_P[COG] = float("{:.2f}".format((num/num_Coding_COGs) * 100))
    for COG, num in Non_Coding_COGs.items():
        Non_Coding_COGs_P[COG] = float("{:.2f}".format((num/num_Non_Coding_COGs) * 100))

    print("Coding -  Y and R Removed: ")
    print(chisquare([Coding_COGs_P['A'], Coding_COGs_P['B'], Coding_COGs_P['C'], Coding_COGs_P['D'], Coding_COGs_P['E'],
                     Coding_COGs_P['F'], Coding_COGs_P['G'], Coding_COGs_P['H'], Coding_COGs_P['I'],
                     Coding_COGs_P['J'], Coding_COGs_P['K'], Coding_COGs_P['L'], Coding_COGs_P['M'], Coding_COGs_P['N'],
                     Coding_COGs_P['O'], Coding_COGs_P['P'], Coding_COGs_P['Q'], Coding_COGs_P['T'],
                     Coding_COGs_P['U'], Coding_COGs_P['Z'],
                     Coding_COGs_P['S']],
                    f_exp=[Gene_COGs_P['A'], Gene_COGs_P['B'], Gene_COGs_P['C'], Gene_COGs_P['D'], Gene_COGs_P['E'],
                           Gene_COGs_P['F'], Gene_COGs_P['G'], Gene_COGs_P['H'], Gene_COGs_P['I'],
                           Gene_COGs_P['J'], Gene_COGs_P['K'], Gene_COGs_P['L'], Gene_COGs_P['M'], Gene_COGs_P['N'],
                           Gene_COGs_P['O'], Gene_COGs_P['P'], Gene_COGs_P['Q'], Gene_COGs_P['T'],
                           Gene_COGs_P['U'], Gene_COGs_P['Z'], Gene_COGs_P['S']]))

######### Removed R and Y


    print(Coding_COGs_P['A'],Coding_COGs_P['B'],Coding_COGs_P['C'],Coding_COGs_P['D'],Coding_COGs_P['E'],Coding_COGs_P['F'],Coding_COGs_P['G'],Coding_COGs_P['H'],Coding_COGs_P['I'],
                     Coding_COGs_P['J'],Coding_COGs_P['K'],Coding_COGs_P['L'],Coding_COGs_P['M'],Coding_COGs_P['N'],Coding_COGs_P['O'],Coding_COGs_P['P'],Coding_COGs_P['Q'],Coding_COGs_P['T'],
                     Coding_COGs_P['U'],Coding_COGs_P['Z'],Coding_COGs_P['S'])

    print(Gene_COGs_P['A'],Gene_COGs_P['B'],Gene_COGs_P['C'],Gene_COGs_P['D'],Gene_COGs_P['E'],Gene_COGs_P['F'],Gene_COGs_P['G'],Gene_COGs_P['H'],Gene_COGs_P['I'],
                     Gene_COGs_P['J'],Gene_COGs_P['K'],Gene_COGs_P['L'],Gene_COGs_P['M'],Gene_COGs_P['N'],Gene_COGs_P['O'],Gene_COGs_P['P'],Gene_COGs_P['Q'],Gene_COGs_P['T'],
                     Gene_COGs_P['U'],Gene_COGs_P['Z'],Gene_COGs_P['S'])

    print("Non Coding - Y,R and B Removed: ") # Removed B for Non Coding
    print(chisquare([Non_Coding_COGs_P['A'], Non_Coding_COGs_P['C'], Non_Coding_COGs_P['D'], Non_Coding_COGs_P['E'],
                     Non_Coding_COGs_P['F'], Non_Coding_COGs_P['G'], Non_Coding_COGs_P['H'], Non_Coding_COGs_P['I'],
                     Non_Coding_COGs_P['J'], Non_Coding_COGs_P['K'], Non_Coding_COGs_P['L'], Non_Coding_COGs_P['M'], Coding_COGs_P['N'],
                     Non_Coding_COGs_P['O'], Non_Coding_COGs_P['P'], Non_Coding_COGs_P['Q'], Non_Coding_COGs_P['T'],
                     Non_Coding_COGs_P['U'], Non_Coding_COGs_P['Z'],
                     Non_Coding_COGs_P['S']],
                    f_exp=[Gene_COGs_P['A'], Gene_COGs_P['C'], Gene_COGs_P['D'], Gene_COGs_P['E'],
                           Gene_COGs_P['F'], Gene_COGs_P['G'], Gene_COGs_P['H'], Gene_COGs_P['I'],
                           Gene_COGs_P['J'], Gene_COGs_P['K'], Gene_COGs_P['L'], Gene_COGs_P['M'], Gene_COGs_P['N'],
                           Gene_COGs_P['O'], Gene_COGs_P['P'], Gene_COGs_P['Q'], Gene_COGs_P['T'],
                           Gene_COGs_P['U'], Gene_COGs_P['Z'], Gene_COGs_P['S']]))

    print(Non_Coding_COGs_P['A'], Non_Coding_COGs_P['C'], Non_Coding_COGs_P['D'], Non_Coding_COGs_P['E'],
          Non_Coding_COGs_P['F'], Non_Coding_COGs_P['G'], Non_Coding_COGs_P['H'], Non_Coding_COGs_P['I'],
          Non_Coding_COGs_P['J'], Non_Coding_COGs_P['K'], Non_Coding_COGs_P['L'], Non_Coding_COGs_P['M'], Non_Coding_COGs_P['N'],
          Non_Coding_COGs_P['O'], Non_Coding_COGs_P['P'], Non_Coding_COGs_P['Q'], Non_Coding_COGs_P['T'],
          Non_Coding_COGs_P['U'], Non_Coding_COGs_P['Z'], Non_Coding_COGs_P['S'])

    print(Gene_COGs_P['A'], Gene_COGs_P['C'], Gene_COGs_P['D'], Gene_COGs_P['E'], Gene_COGs_P['F'],
          Gene_COGs_P['G'], Gene_COGs_P['H'], Gene_COGs_P['I'],
          Gene_COGs_P['J'], Gene_COGs_P['K'], Gene_COGs_P['L'], Gene_COGs_P['M'], Gene_COGs_P['N'], Gene_COGs_P['O'],
          Gene_COGs_P['P'], Gene_COGs_P['Q'], Gene_COGs_P['T'],
          Gene_COGs_P['U'], Gene_COGs_P['Z'], Gene_COGs_P['S'])




