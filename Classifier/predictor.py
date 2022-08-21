import collections
from utils import *
import sys
import gc


def predictor(Reads,model,options):
    multi_Coding_Reads = collections.defaultdict(list)

    params = get_params(11)
    params["pool_size"] = 1
    MAXLEN = 75
    batch_size = 50
    id2seq = {}
    aaletters = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'}
    aaindex = {'K': 0, 'H': 1, 'T': 2, 'V': 3, 'Q': 4, 'L': 5, 'D': 6, 'R': 7, 'A': 8, 'G': 9, 'Y': 10, 'F': 11,
               'S': 12, 'N': 13, 'P': 14, 'C': 15, 'M': 16, 'I': 17, 'W': 18, 'E': 19}

    for Read, Frames in Reads.items():
        for Frame in Frames:
            id2seq[Frame[0]] = Frame[1]

    prot2embed = {}
    for idx in id2seq:
        prot2embed[idx] = to_onehot(id2seq[idx], aaindex)

    tuple_run = {x: 1 for x in id2seq} # Convert data into format needed for model
    tuple_run = list(tuple_run.items())
    tuple_run = np.array(tuple_run)
    data_generator = seq_Generator(tuple_run[:,0], tuple_run[:,1], batch_size, prot2embed, MAXLEN=MAXLEN)

    classifying = model.predict_generator(data_generator)

    idfromseq = list(id2seq.keys())
    classif = classifying.tolist()
    classif = [float(i[0]) for i in classif]

    classified = list(zip(idfromseq, classif))

    ## Turn classified into a dict or loop through it and turn it into a dict by picking most likely coding frame by using the frame with highest score
    Coding = open(options.out_prefix+'_Coding.fa','a')
    Non_Coding = open(options.out_prefix+'_Non_Coding.fa', 'a')
    for key, value in classified:
        current_read = key.split('_')[0]
        frames = Reads[current_read]
        for i, frame in enumerate(frames):
            if key in frame:
                Reads[current_read][i].append(value)
                break

    for Read, Frames in Reads.items():
        highest_scored_frame = max(Frames, key=lambda x: x[2]) # Frame with highest score
        if highest_scored_frame[2] >= 0.5:
            Coding.write(highest_scored_frame[0]+'_Score:'+ str(format(highest_scored_frame[2], '.2f')) + '\n' +highest_scored_frame[1]+'\n')
        for Frame in Frames:
        ##### Need to add stops back in
            if Frame[2] < 0.5:
                Non_Coding.write(Frame[0] + '_Score:' + str(format(Frame[2], '.2f')) + '\n' + Frame[1] + '\n')
            elif Frame == highest_scored_frame:
                continue
            elif Frame != highest_scored_frame and Frame[2] >= 0.5:
                multi_Coding_Reads[Read].append([highest_scored_frame, Frame])
                Coding.write(Frame[0] + '_Score:' + str(format(Frame[2], '.2f')) + '\n' + Frame[1] + '\n')
            elif Frame != highest_scored_frame and Frame[2] < 0.5:
                Non_Coding.write(Frame[0] +'_Score:'+ str(format(Frame[2], '.2f')) + '\n' + Frame[1] + '\n')
            else:
                sys.exit("Crashed")

    # multi_out = open('./multi_out.fa','a')
    # for key, value in multi_Coding_Reads.items():
    #     multi_out.write(str(key)+'_'+str(value)+'\n')

    del classified, classifying, Reads, idfromseq, id2seq, model, multi_Coding_Reads, tuple_run, data_generator, prot2embed
    gc.collect()
    print("Done")
