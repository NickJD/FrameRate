import argparse
import collections
import numpy as np
import sys


def remove_frame_with_None(collection_list):
    none_filtered_collection_list = []
    for list in collection_list:
        if None not in list:
            none_filtered_collection_list.append(list)

    return none_filtered_collection_list

def calculate_parameters(collection_lengths,collection_scores,collection_pident,collection_bitscores,):
    ###################################
    num_frames = len(collection_lengths)
    mean_length = round(np.mean(collection_lengths),2)
    percentile_length_25 = np.percentile(collection_lengths,25)
    percentile_length_75 = np.percentile(collection_lengths, 75)
    std_length = round(np.std(collection_lengths),2)
    min_length = np.min(collection_lengths)
    max_length = np.max(collection_lengths)
    ####################################
    mean_score = round(np.mean(collection_scores),2)
    percentile_score_25 = np.percentile(collection_scores,25)
    percentile_score_75 = np.percentile(collection_scores, 75)
    std_score = round(np.std(collection_scores),2)
    min_score = np.min(collection_scores)
    max_score = np.max(collection_scores)
    ####################################
    if collection_pident == None and collection_bitscores == None:
        return num_frames,mean_length,percentile_length_25,percentile_length_75,std_length,min_length,max_length,\
               mean_score,percentile_score_25,percentile_score_75,std_score,min_score,max_score
    else:
        ####################################
        mean_pident = round(np.mean(collection_pident), 2)
        percentile_pident_25 = np.percentile(collection_pident, 25)
        percentile_pident_75 = np.percentile(collection_pident, 75)
        std_pident = round(np.std(collection_pident), 2)
        min_pident = np.min(collection_pident)
        max_pident = np.max(collection_pident)
        ####################################
        mean_bitscores = round(np.mean(collection_bitscores), 2)
        percentile_bitscores_25 = np.percentile(collection_bitscores, 25)
        percentile_bitscores_75 = np.percentile(collection_bitscores, 75)
        std_bitscores = round(np.std(collection_bitscores), 2)
        min_bitscores = np.min(collection_bitscores)
        max_bitscores = np.max(collection_bitscores)
        ####################################
        return num_frames,mean_length,percentile_length_25,percentile_length_75,std_length,min_length,max_length,\
               mean_score,percentile_score_25,percentile_score_75,std_score,min_score,max_score,mean_pident,percentile_pident_25,\
               percentile_pident_75,std_pident,min_pident,max_pident,mean_bitscores,percentile_bitscores_25,percentile_bitscores_75,\
               std_bitscores,min_bitscores,max_bitscores




def frame_scores(frames,collection):
    current_id = ''
    for line in frames:
        if line.startswith('>'):
            score = line.split('Score:')[1].strip()
            collection[line.split(' ')[0].strip().strip('>')].append(float(score)) # reove the split here
            current_id = line.split(' ')[0].strip().strip('>') # and here
        else:
            collection[current_id].append(len(line.strip()))

    collection_list = list(collection.values())
    collection_scores = [el[0] for el in collection_list]
    collection_lengths = [el[1] for el in collection_list]

    num_frames, mean_length, percentile_length_25, percentile_length_75, std_length, min_length, max_length, \
    mean_score, percentile_score_25, percentile_score_75, std_score, min_score, max_score = calculate_parameters(collection_lengths,collection_scores,None,None)

    return collection,num_frames,mean_length,percentile_length_25,percentile_length_75,std_length,min_length,max_length,\
           mean_score,percentile_score_25,percentile_score_75,std_score,min_score,max_score

def swiss_scores(frames,collection):
    for line in frames:
        line_data = line.split('\t')
        collection[line_data[0]].append(float(line_data[2])) # PIdent
        collection[line_data[0]].append(float(line_data[11].strip('\n'))) # BitScore
    for frame, data in collection.items():
        if len(data) == 2:
            collection[frame].extend([None,None])  # Adding None to those frames without a hit

    collection_list = list(collection.values())
    none_filtered_collection_list = remove_frame_with_None(collection_list)
    collection_scores = [el[0] for el in none_filtered_collection_list]
    collection_lengths = [el[1] for el in none_filtered_collection_list]
    collection_pident = [el[2] for el in none_filtered_collection_list]
    collection_bitscore = [el[3] for el in none_filtered_collection_list]
    num_frames, mean_length, percentile_length_25, percentile_length_75, std_length, min_length, max_length, \
    mean_score, percentile_score_25, percentile_score_75, std_score, min_score, max_score, mean_pident, percentile_pident_25, \
    percentile_pident_75, std_pident, min_pident, max_pident, mean_bitscores, percentile_bitscores_25, percentile_bitscores_75, \
    std_bitscores, min_bitscores, max_bitscores = calculate_parameters(collection_lengths, collection_scores,collection_pident,collection_bitscore)


    return collection,num_frames,mean_length,percentile_length_25,percentile_length_75,std_length,min_length,max_length,\
               mean_score,percentile_score_25,percentile_score_75,std_score,min_score,max_score,mean_pident,percentile_pident_25,\
               percentile_pident_75,std_pident,min_pident,max_pident,mean_bitscores,percentile_bitscores_25,percentile_bitscores_75,\
               std_bitscores,min_bitscores,max_bitscores


def EggNOG_frame_scores(frames,collection,COGs):
    num_hits = [0, 0]
    frames_with_hit_count = 0
    for line in frames:
        if line.startswith('#'):
            continue
        else:
            line = line.split('\t')
            COG = line[6]
            COGs[COG]+=1
            frames_with_hit_count +=1
            collection[line[0]].extend(COG)

    num_hits[0] = frames_with_hit_count
    num_hits[1] = COGs['-']

    for frame, data in collection.items():
        if len(data) == 4:
            collection[frame].extend([None])  # Adding None to those frames without a hit
    collection_list = list(collection.values())
    none_filtered_collection_list = remove_frame_with_None(collection_list)
    collection_scores = [el[0] for el in none_filtered_collection_list]
    collection_lengths = [el[1] for el in none_filtered_collection_list]
    collection_pident = [el[2] for el in none_filtered_collection_list]
    collection_bitscore = [el[3] for el in none_filtered_collection_list]
    num_frames, mean_length, percentile_length_25, percentile_length_75, std_length, min_length, max_length, \
    mean_score, percentile_score_25, percentile_score_75, std_score, min_score, max_score, mean_pident, percentile_pident_25, \
    percentile_pident_75, std_pident, min_pident, max_pident, mean_bitscores, percentile_bitscores_25, percentile_bitscores_75, \
    std_bitscores, min_bitscores, max_bitscores = calculate_parameters(collection_lengths, collection_scores,collection_pident,collection_bitscore)

    return collection, COGs,num_hits, num_frames,mean_length,percentile_length_25,percentile_length_75,std_length,min_length,max_length,\
               mean_score,percentile_score_25,percentile_score_75,std_score,min_score,max_score,mean_pident,percentile_pident_25,\
               percentile_pident_75,std_pident,min_pident,max_pident,mean_bitscores,percentile_bitscores_25,percentile_bitscores_75,\
               std_bitscores,min_bitscores,max_bitscores

def EggNOG_grouping(COGs):
    COG_groups = {'Cell':0,'Info':0,'Meta':0,'Poor':0}
    Cell = ['D','M','N','O','T','U','V','W','Y','Z']
    Info = ['A','B','J','K','L']
    Meta = ['C','E','F','G','H','I','P','Q']
    Poor = ['R','S']

    for COG, count in COGs.items():
        if COG in Cell:
            COG_groups['Cell'] +=count
        elif COG in Info:
            COG_groups['Info'] +=count
        elif COG in Meta:
            COG_groups['Meta'] +=count
        elif COG in Poor:
            COG_groups['Poor'] +=count
        #else:
            #print("Unknown COG: " + str(COG))
    return COG_groups


def EggNOG_CDS_scores(genes,collection,COGs):
    num_hits = [0,0]
    frames_with_hit_count = 0
    for line in genes:
        if line.startswith('#'):
            continue
        else:
            line = line.split('\t')
            COG = line[6]
            COGs[COG]+=1
            frames_with_hit_count +=1
            collection[line[0]].extend(COG)
    num_hits[0] = frames_with_hit_count
    num_hits[1] = COGs['-']


    return collection,COGs,num_hits

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', action='store', dest='coding', required=True,
                        help='Coding Frames')
    parser.add_argument('-nc', action='store', dest='non_coding', required=False,
                        help='Non_Coding Frames')
    parser.add_argument('-swiss_c', action='store', dest='swiss_coding', required=False,
                        help='SwissProt Coding Frames')
    parser.add_argument('-swiss_nc', action='store', dest='swiss_non_coding', required=False,
                        help='SwissProt Non_Coding Frames')
    parser.add_argument('-eggnog_c', action='store', dest='eggnog_coding', required=False,
                        help='Eggnog Coding Frames')
    parser.add_argument('-eggnog_nc', action='store', dest='eggnog_non_coding', required=False,
                        help='Eggnog Non_Coding Frames')
    parser.add_argument('-cds', action='store', dest='eggnog_cds', required=False,
                        help='Eggnog CDS Genes')

    options = parser.parse_args()

    Coding_Frames_Collection = collections.defaultdict(list)
    Non_Coding_Frame_Collection = collections.defaultdict(list)
    ################## Reads
    Coding_Frames = open(options.coding, 'r')
    Coding_Frames_Collection,c_num_frames,c_mean_length,c_percentile_length_25,c_percentile_length_75,c_std_length,c_min_length,c_max_length,\
           c_mean_score,c_percentile_score_25,c_percentile_score_75,c_std_score,c_min_score,c_max_score = frame_scores(Coding_Frames,Coding_Frames_Collection)
    ###
    Non_Coding_Frames = open(options.non_coding, 'r')
    Non_Coding_Frames_Collection,nc_num_frames,nc_mean_length,nc_percentile_length_25,nc_percentile_length_75,nc_std_length,nc_min_length,nc_max_length,\
           nc_mean_score,nc_percentile_score_25,nc_percentile_score_75,nc_std_score,nc_min_score,nc_max_score = frame_scores(Non_Coding_Frames,Non_Coding_Frame_Collection)
    print("FrameRate Output Frame Metrics:")
    print("Metric\tCoding Frames\tNon-Coding Frames")
    print("Number of Frames\t"+str(c_num_frames)+"\t"+str(nc_num_frames)+"\nMean confidence score\t"+str(c_mean_score)+"\t"+str(nc_mean_score)+"\n75th-Percentile score\t"+
          str(c_percentile_score_75)+"\t"+str(nc_percentile_score_75)+"\n25th-Percentile score\t"+str(c_percentile_score_25)+"\t"+str(nc_percentile_score_25)+"\nStandard Dev score\t"+
          str(c_std_score)+"\t"+str(nc_std_score)+"\nMinimum confidence score\t"+str(c_min_score)+"\t"+str(nc_min_score)+"\nMaximum confidence score\t"+str(c_max_score)+"\t"+str(nc_max_score)+
          "\nMean length\t"+str(c_mean_length)+"\t"+str(nc_mean_length)+"\n75th-Percentile length\t"+str(c_percentile_length_75)+"\t"+str(nc_percentile_length_75)+"\n25th-Percentile length\t"+
          str(c_percentile_length_25)+"\t"+str(nc_percentile_length_25)+"\nStandard Dev length\t"+str(c_std_length)+"\t"+str(nc_std_length)+"\nMinimum length\t"+str(c_min_length)+"\t"+str(nc_min_length)+
          "\nMaximum length\t"+str(c_max_length)+"\t"+str(nc_max_length))
    ################### Swiss-Prot

    Coding_Frames_Swiss = open(options.swiss_coding, 'r')
    Coding_Frames_Collection,c_swiss_num_frames,c_swiss_mean_length,c_swiss_percentile_length_25,c_swiss_percentile_length_75,c_swiss_std_length,\
    c_swiss_min_length,c_swiss_max_length,c_swiss_mean_score,c_swiss_percentile_score_25,c_swiss_percentile_score_75,c_swiss_std_score,c_swiss_min_score,\
    c_swiss_max_score,c_swiss_mean_pident,c_swiss_percentile_pident_25,c_swiss_percentile_pident_75,c_swiss_std_pident,c_swiss_min_pident,c_swiss_max_pident,\
    c_swiss_mean_bitscores,c_swiss_percentile_bitscores_25,c_swiss_percentile_bitscores_75,c_swiss_std_bitscores,c_swiss_min_bitscores,c_swiss_max_bitscores \
        = swiss_scores(Coding_Frames_Swiss,Coding_Frames_Collection)
    ###
    Non_Coding_Frames_Swiss = open(options.swiss_non_coding, 'r')
    Non_Coding_Frame_Collection,nc_swiss_num_frames,nc_swiss_mean_length,nc_swiss_percentile_length_25,nc_swiss_percentile_length_75,nc_swiss_std_length,\
    nc_swiss_min_length,nc_swiss_max_length,nc_swiss_mean_score,nc_swiss_percentile_score_25,nc_swiss_percentile_score_75,nc_swiss_std_score,nc_swiss_min_score,\
    nc_swiss_max_score,nc_swiss_mean_pident,nc_swiss_percentile_pident_25,nc_swiss_percentile_pident_75,nc_swiss_std_pident,nc_swiss_min_pident,nc_swiss_max_pident,\
    nc_swiss_mean_bitscores,nc_swiss_percentile_bitscores_25,nc_swiss_percentile_bitscores_75,nc_swiss_std_bitscores,nc_swiss_min_bitscores,nc_swiss_max_bitscores \
        = swiss_scores(Non_Coding_Frames_Swiss,Non_Coding_Frame_Collection)
    print("\n\nFrameRate Output Frame with Swiss-Prot Hit Metrics:")
    print("Metric\tCoding Frames\tNon-Coding Frames")
    print("Number of Frames\t"+str(c_swiss_num_frames)+"\t"+str(nc_swiss_num_frames)+"\nMean confidence score\t"+str(c_swiss_mean_score)+"\t"+str(nc_swiss_mean_score)+"\n75th-Percentile score\t"+
          str(c_swiss_percentile_score_75)+"\t"+str(nc_swiss_percentile_score_75)+"\n25th-Percentile score\t"+str(c_swiss_percentile_score_25)+"\t"+str(nc_swiss_percentile_score_25)+"\nStandard Dev score\t"+
          str(c_swiss_std_score)+"\t"+str(nc_swiss_std_score)+"\nMinimum confidence score\t"+str(c_swiss_min_score)+"\t"+str(nc_swiss_min_score)+"\nMaximum confidence score\t"+str(c_swiss_max_score)+"\t"+
          str(nc_swiss_max_score)+"\nMean length\t"+str(c_swiss_mean_length)+"\t"+str(nc_swiss_mean_length)+"\n75th-Percentile score\t"+str(c_swiss_percentile_length_75)+"\t"+str(nc_swiss_percentile_length_75)+
          "\n25th-Percentile score\t"+str(c_swiss_percentile_length_25)+"\t"+str(nc_swiss_percentile_length_25)+"\nStandard Dev length\t"+str(c_swiss_std_length)+"\t"+str(nc_swiss_std_length)+"\nMinimum length\t"+
          str(c_swiss_min_length)+"\t"+str(nc_swiss_min_length)+"\nMaximum length\t"+str(c_swiss_max_length)+"\t"+str(nc_swiss_max_length)+

          "\nMean PIdent score\t"+str(c_swiss_mean_pident)+"\t"+str(nc_swiss_mean_pident)+"\n75th-Percentile PIdent score\t"+str(c_swiss_percentile_pident_75)+"\t"+str(nc_swiss_percentile_pident_75)+
          "\n25th-Percentile PIdent score\t"+str(c_swiss_percentile_pident_25)+"\t"+str(nc_swiss_percentile_pident_25)+"\nStandard Dev PIdent score\t"+str(c_swiss_std_pident)+"\t"+str(nc_swiss_std_pident)+
          "\nMinimum confidence PIdent score\t"+str(c_swiss_min_pident)+"\t"+str(nc_swiss_min_pident)+"\nMaximum confidence PIdent score\t"+str(c_swiss_max_pident)+"\t"+

          "\nMean Bitscore\t"+str(c_swiss_mean_bitscores)+"\t"+str(nc_swiss_mean_bitscores)+"\n75th-Percentile Bitscore\t"+str(c_swiss_percentile_bitscores_75)+"\t"+str(nc_swiss_percentile_bitscores_75)+
          "\n25th-Percentile Bitscore\t"+str(c_swiss_percentile_bitscores_25)+"\t"+str(nc_swiss_percentile_bitscores_25)+"\nStandard Dev Bitscore\t"+str(c_swiss_std_bitscores)+"\t"+str(nc_swiss_std_bitscores)+
          "\nMinimum confidence Bitscore\t"+str(c_swiss_min_bitscores)+"\t"+str(nc_swiss_min_bitscores)+"\nMaximum confidence Bitscore\t"+str(c_swiss_max_bitscores)+"\t"+str(nc_swiss_max_bitscores))
    ################### EggNOG
    print("\nEggNOG")
    Coding_Frames_EggNOG = open(options.eggnog_coding, 'r')
    Coding_COGs = collections.defaultdict(int)
    Coding_Frames_Collection, Coding_COGs, c_eggnog_num_hits,c_eggnog_num_frames,c_eggnog_mean_length,c_eggnog_percentile_length_25,c_eggnog_percentile_length_75,c_eggnog_std_length,\
    c_eggnog_min_length,c_eggnog_max_length,c_eggnog_mean_score,c_eggnog_percentile_score_25,c_eggnog_percentile_score_75,c_eggnog_std_score,c_eggnog_min_score,c_eggnog_max_score,c_eggnog_mean_pident,\
    c_eggnog_percentile_pident_25,c_eggnog_percentile_pident_75,c_eggnog_std_pident,c_eggnog_min_pident,c_eggnog_max_pident,c_eggnog_mean_bitscores,c_eggnog_percentile_bitscores_25,c_eggnog_percentile_bitscores_75,\
    c_eggnog_std_bitscores,c_eggnog_min_bitscores,c_eggnog_max_bitscores = EggNOG_frame_scores(Coding_Frames_EggNOG, Coding_Frames_Collection,Coding_COGs)
    ###
    Non_Coding_Frames_EggNOG = open(options.eggnog_non_coding, 'r')
    Non_Coding_COGs = collections.defaultdict(int)
    Non_Coding_Frames_Collection, Non_Coding_COGs, nc_eggnog_num_hits,nc_eggnog_num_frames,nc_eggnog_mean_length,nc_eggnog_percentile_length_25,nc_eggnog_percentile_length_75,nc_eggnog_std_length,\
    nc_eggnog_min_length,nc_eggnog_max_length,nc_eggnog_mean_score,nc_eggnog_percentile_score_25,nc_eggnog_percentile_score_75,nc_eggnog_std_score,nc_eggnog_min_score,nc_eggnog_max_score,nc_eggnog_mean_pident,\
    nc_eggnog_percentile_pident_25,nc_eggnog_percentile_pident_75,nc_eggnog_std_pident,nc_eggnog_min_pident,nc_eggnog_max_pident,nc_eggnog_mean_bitscores,nc_eggnog_percentile_bitscores_25,nc_eggnog_percentile_bitscores_75,\
    nc_eggnog_std_bitscores,nc_eggnog_min_bitscores,nc_eggnog_max_bitscores = EggNOG_frame_scores(Non_Coding_Frames_EggNOG, Non_Coding_Frames_Collection,Non_Coding_COGs)
    ###
    print("FrameRate Output Frame with EggNOG Hit Metrics:")
    print("Metric\tCoding Frames\tNon-Coding Frames")
    print("Number of Frames\t"+str(c_eggnog_num_frames)+"\t"+str(nc_eggnog_num_frames)+"\nMean confidence score\t"+str(c_eggnog_mean_score)+"\t"+str(nc_eggnog_mean_score)+"\n75th-Percentile score\t"+
          str(c_eggnog_percentile_score_75)+"\t"+str(nc_eggnog_percentile_score_75)+"\n25th-Percentile score\t"+str(c_eggnog_percentile_score_25)+"\t"+str(nc_eggnog_percentile_score_25)+"\nStandard Dev score\t"+
          str(c_eggnog_std_score)+"\t"+str(c_eggnog_std_score)+"\nMinimum confidence score\t"+str(c_eggnog_min_score)+"\t"+str(nc_eggnog_min_score)+"\nMaximum confidence score\t"+str(c_eggnog_max_score)+"\t"+
          str(nc_eggnog_max_score)+"\nMean length\t"+str(c_eggnog_mean_length)+"\t"+str(nc_eggnog_mean_length)+"\n75th-Percentile score\t"+str(c_eggnog_percentile_length_75)+"\t"+str(nc_eggnog_percentile_length_75)+
          "\n25th-Percentile score\t"+str(c_eggnog_percentile_length_25)+"\t"+str(nc_eggnog_percentile_length_25)+"\nStandard Dev length\t"+str(c_eggnog_std_length)+"\t"+str(nc_eggnog_std_length)+"\nMinimum length\t"+
          str(c_eggnog_min_length)+"\t"+str(nc_eggnog_min_length)+"\nMaximum length\t"+str(c_eggnog_max_length)+"\t"+str(nc_eggnog_max_length))

    #####
    Coding_COG_groups = EggNOG_grouping(Coding_COGs)
    Non_Coding_COG_groups = EggNOG_grouping(Non_Coding_COGs)
    ###################
    EggNOG_CDS = open(options.eggnog_cds, 'r')
    CDS_COGs = collections.defaultdict(int)
    EggNOG_CDS_Collection = collections.defaultdict(list)
    EggNOG_CDS_Collection,CDS_COGs,CDS_num_hits = EggNOG_CDS_scores(EggNOG_CDS, EggNOG_CDS_Collection,CDS_COGs)
    CDS_COG_groups = EggNOG_grouping(CDS_COGs)
    ##################
    Coding_COG_group_Cell_p = (Coding_COG_groups['Cell']/sum(Coding_COG_groups.values())) * 100
    Coding_COG_group_Info_p = (Coding_COG_groups['Info'] / sum(Coding_COG_groups.values())) * 100
    Coding_COG_group_Meta_p = (Coding_COG_groups['Meta'] / sum(Coding_COG_groups.values())) * 100
    Coding_COG_group_Poor_p = (Coding_COG_groups['Poor'] / sum(Coding_COG_groups.values())) * 100
    ###
    Non_Coding_COG_group_Cell_p = (Non_Coding_COG_groups['Cell']/sum(Non_Coding_COG_groups.values())) * 100
    Non_Coding_COG_group_Info_p = (Non_Coding_COG_groups['Info'] / sum(Non_Coding_COG_groups.values())) * 100
    Non_Coding_COG_group_Meta_p = (Non_Coding_COG_groups['Meta'] / sum(Non_Coding_COG_groups.values())) * 100
    Non_Coding_COG_group_Poor_p = (Non_Coding_COG_groups['Poor'] / sum(Non_Coding_COG_groups.values())) * 100
    ###
    CDS_COG_group_Cell_p = (CDS_COG_groups['Cell']/sum(CDS_COG_groups.values())) * 100
    CDS_COG_group_Info_p = (CDS_COG_groups['Info'] / sum(CDS_COG_groups.values())) * 100
    CDS_COG_group_Meta_p = (CDS_COG_groups['Meta'] / sum(CDS_COG_groups.values())) * 100
    CDS_COG_group_Poor_p = (CDS_COG_groups['Poor'] / sum(CDS_COG_groups.values())) * 100
    ###################
    print("\nCOGs:")
    print("Metric\tCoding Frames\tNon-Coding Frames\tCDS Genes")
    try:
        print("Number of sequences with an EggNOG Hit | With a COG | %\t"+str(c_eggnog_num_hits[0])+"|"+str(c_eggnog_num_hits[1])+"|"+str(round((c_eggnog_num_hits[1]/c_eggnog_num_frames)*100,2))+"\t"+
        str(nc_eggnog_num_hits[1])+"|"+str(nc_eggnog_num_hits[0])+"|"+str(round((nc_eggnog_num_hits[1]/nc_eggnog_num_frames)*100,2))+"\t"+
        str(CDS_num_hits[1])+"|"+str(CDS_num_hits[0])+"|"+str(round((CDS_num_hits[1]/CDS_num_hits)*100,2)))
        print("CELLULAR PROCESSES AND SIGNALING:\t" + str(Coding_COG_groups['Cell'])+"|"+str(round(Coding_COG_group_Cell_p,2))+"\t"+str(Non_Coding_COG_groups['Cell'])+"|"+str(round(Non_Coding_COG_group_Cell_p,2))+"\t"+str(CDS_COG_groups['Cell'])+"|"+str(round(CDS_COG_group_Cell_p,2)))
        print("D\t"+str(Coding_COGs['D'])+"\t"+str(Non_Coding_COGs['D'])+"\t"+str(CDS_COGs['D'])+"\nM\t"+str(Coding_COGs['M'])+"\t"+str(Non_Coding_COGs['M'])+
              "\t"+str(CDS_COGs['M'])+"\nN\t"+str(Coding_COGs['N'])+"\t"+str(Non_Coding_COGs['N'])+"\t"+str(CDS_COGs['N'])+"\nO\t"+str(Coding_COGs['O'])+"\t"+
              str(Non_Coding_COGs['O'])+"\t"+str(CDS_COGs['O'])+"\nT\t"+str(Coding_COGs['T'])+"\t"+str(Non_Coding_COGs['T'])+"\t"+str(CDS_COGs['T'])+"\nU\t"+
              str(Coding_COGs['U'])+"\t"+str(Coding_COGs['U'])+"\t"+str(CDS_COGs['U'])+"\nV\t"+str(Coding_COGs['V'])+"\t"+str(Coding_COGs['V'])+"\t"+str(CDS_COGs['V'])+
              "\nW\t"+str(Coding_COGs['W'])+"\t"+str(Coding_COGs['W'])+"\t"+str(CDS_COGs['W'])+"\nY\t"+str(Coding_COGs['Y'])+"\t"+str(Coding_COGs['Y'])+"\t"+str(CDS_COGs['Y'])+
              "\nZ\t"+str(Coding_COGs['Z'])+"\t"+str(Coding_COGs['Z'])+"\t"+str(CDS_COGs['Z']))
        print("INFORMATION STORAGE AND PROCESSING\t" + str(Coding_COG_groups['Info'])+"|"+str(round(Coding_COG_group_Info_p,2))+"\t"+str(Non_Coding_COG_groups['Info'])+"|"+str(round(Non_Coding_COG_group_Info_p,2))+"\t"+str(CDS_COG_groups['Info'])+"|"+str(round(CDS_COG_group_Info_p,2)))
        print("A\t"+str(Coding_COGs['A'])+"\t"+str(Non_Coding_COGs['A'])+"\t"+str(CDS_COGs['A'])+"\nB\t"+str(Coding_COGs['B'])+"\t"+str(Non_Coding_COGs['B'])+
              "\t"+str(CDS_COGs['B'])+"\nJ\t"+str(Coding_COGs['J'])+"\t"+str(Non_Coding_COGs['J'])+"\t"+str(CDS_COGs['J'])+"\nK\t"+str(Coding_COGs['K'])+"\t"+
              str(Non_Coding_COGs['K'])+"\t"+str(CDS_COGs['K'])+"\nL\t" + str(Coding_COGs['L']) + "\t" + str(Non_Coding_COGs['L'])+"\t"+str(CDS_COGs['L']))
        print("METABOLISM\t" + str(Coding_COG_groups['Meta'])+"|"+str(round(Coding_COG_group_Meta_p,2))+"\t"+str(Non_Coding_COG_groups['Meta'])+"|"+str(round(Non_Coding_COG_group_Meta_p,2))+"\t"+str(CDS_COG_groups['Meta'])+"|"+str(round(CDS_COG_group_Meta_p,2)))
        print("C\t"+str(Coding_COGs['C'])+"\t"+str(Non_Coding_COGs['C'])+"\t"+str(CDS_COGs['C'])+"\nE\t"+str(Coding_COGs['E'])+"\t"+str(Non_Coding_COGs['E'])+
              "\t"+str(CDS_COGs['E'])+"\nF\t"+str(Coding_COGs['F'])+"\t"+str(Non_Coding_COGs['F'])+"\t"+str(CDS_COGs['F'])+"\nG\t"+str(Coding_COGs['G'])+
              "\t"+str(Non_Coding_COGs['G'])+"\t"+str(CDS_COGs['G'])+"\nH\t" + str(Coding_COGs['H']) + "\t" + str(Non_Coding_COGs['H'])+"\t"+str(CDS_COGs['H'])+
              "\nI\t"+str(Coding_COGs['I'])+"\t"+str(Non_Coding_COGs['I'])+"\t"+str(CDS_COGs['I'])+"\nP\t"+str(Coding_COGs['P'])+"\t"+str(Non_Coding_COGs['P'])+
              "\t"+str(CDS_COGs['P'])+"\nQ\t"+str(Coding_COGs['Q'])+"\t"+str(Non_Coding_COGs['Q'])+"\t"+str(CDS_COGs['Q']))
        print("POORLY CHARACTERIZED\t" + str(Coding_COG_groups['Poor'])+"|"+str(round(Coding_COG_group_Poor_p,2))+"\t"+str(Non_Coding_COG_groups['Poor'])+"|"+str(round(Non_Coding_COG_group_Poor_p,2))+"\t"+str(CDS_COG_groups['Poor'])+"|"+str(round(CDS_COG_group_Poor_p,2)))
        print("\nR\t"+str(Coding_COGs['R'])+"\t"+str(Non_Coding_COGs['R'])+"\t"+str(CDS_COGs['R'])+"\nS\t"+str(Coding_COGs['S'])+"\t"+str(Non_Coding_COGs['S'])+"\t"+str(CDS_COGs['S']))
    except ZeroDivisionError:
        sys.exit("ZeroDivisionError: division by zero")