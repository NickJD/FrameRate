import numpy as np
from tensorflow.python.client import device_lib
import tensorflow.compat.v1 as tf
import random
from utils import *
from random import seed
from random import randint
import sys
params = get_params(11)
params["pool_size"] = 1
#params["max_kernel"] = 17
params["nb_filters"] = 4

print("Params: "+str(params))

random.seed(1) # For random TN selector

MAXLEN = 75
batch_size = 200
window = 75
overlap = 25
positives, negatives = [], []
id2seq = {}
aaletters = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
aaindex = {'K': 0, 'H': 1, 'T': 2, 'V': 3, 'Q': 4, 'L': 5, 'D': 6, 'R': 7, 'A': 8, 'G': 9, 'Y': 10, 'F': 11,
           'S': 12, 'N': 13, 'P': 14, 'C': 15, 'M': 16, 'I': 17, 'W': 18, 'E': 19}
# Test - Move 1-20 to 0-19


model_file = "best_model_Taxa"


#tf.config.set_visible_devices([], 'GPU') # TO Run CPU

print("Setup:")
print(device_lib.list_local_devices())
print("Num GPUs Available: ", len(tf.config.experimental.list_physical_devices('GPU')))



def get_overlapped_chunks(textin, chunksize, overlapsize):  # Maybe remove last chunk as it will be covered and will be short???
    return [ textin[a:a+chunksize] for a in range(0,len(textin), chunksize-overlapsize)]

### check the overlap


################################
with open('/home/nick/Nextcloud/Dimsum/Extended_CoDing_Sequences_For_Training_Biggest.csv', 'r') as f: # Training Data
    count = 0
    for line in f:
        items = line.strip().split(',')
        #if count >= 1000:
        #    break
        if aaletters.issuperset(items[3]): # Actinomycetota
            if 'Mycobacterium' in items[0] and int(items[2]) == 1:
                positives.append((items[0],1))
                id2seq[items[0]] = items[3]
            # elif 'Actinomyces' in items[0] and int(items[2]) == 1:
            #     positives.append((items[0],1))
            #     id2seq[items[0]] = items[3]
            # elif 'Streptomyces' in items[0] and int(items[2]) == 1:
            #     positives.append((items[0],1))
            #     id2seq[items[0]] = items[3]
            # elif 'Gardnerella' in items[0] and int(items[2]) == 1:
            #     positives.append((items[0],1))
            #     id2seq[items[0]] = items[3]
            # elif 'Bifidobacterium' in items[0] and int(items[2]) == 1:
            #     positives.append((items[0], 1))
            #     id2seq[items[0]] = items[3]
            ############################################  # Pseudomonadota / Proteobacteria
            # elif 'Escherichia' in items[0] and int(items[2]) == 1: # randomly - using seed - pick one of 5 negatives for each positive
            #     negatives.append((items[0], 0))
            #     #positives.append((items[0], 1)) # mess with training data
            #     id2seq[items[0]] = items[3]
            # elif 'Shigella' in items[0] and int(items[2]) == 1:
            #     negatives.append((items[0], 0))
            #     id2seq[items[0]] = items[3]
            # elif 'Salmonella' in items[0] and int(items[2]) == 1:
            #     negatives.append((items[0], 0))
            #     id2seq[items[0]] = items[3]
            # elif 'Vibrio' in items[0] and int(items[2]) == 1:
            #     negatives.append((items[0], 0))
            #     id2seq[items[0]] = items[3]
            # elif 'Yersinia' in items[0] and int(items[2]) == 1:
            #     negatives.append((items[0], 0))
            #     id2seq[items[0]] = items[3]
            elif 'Legionella' in items[0] and int(items[2]) == 1:
                 negatives.append((items[0], 0))
                 id2seq[items[0]] = items[3]
            ################ testing!!!



print("Positive-Negative-IDs-AA")
print(len(positives), len(negatives), len(id2seq), len(aaletters))
random.shuffle(positives)
#sub_positives = random.sample(positives, 40000)
random.shuffle(negatives)# Extra shuffling
#sub_negatives = random.sample(negatives, 40000)
#######################################
#np.save('id2seq.npy',id2seq)
#np.save('sub_positives.npy',sub_positives)
#np.save('sub_negatives.npy',sub_negatives)
###################################
# print("Loading")
# sub_positives = np.load('sub_positives.npy')
# sub_negatives = np.load('sub_negatives.npy')
# id2 = np.load('id2seq.npy',allow_pickle=True)
# i = dict(np.ndenumerate(id2))
# id2seq = i[()]
# sub_positives = sub_positives.tolist()
# sub_negatives = sub_negatives.tolist()
#################################
print("Splitting")
train_pos, val_pos, test_pos = split_train(positives, 0.8, 0.9)
train_neg, val_neg, test_neg = split_train(negatives, 0.8, 0.9)
print(len(train_pos), len(val_pos), len(test_pos))
print(len(train_neg), len(val_neg), len(test_neg))

tuple_train = np.concatenate((train_pos, train_neg[:len(train_pos)]), axis=0)
tuple_val = np.concatenate((val_pos, val_neg[:len(val_pos)]), axis = 0)
tuple_test = np.concatenate((test_pos, test_neg[:len(test_pos)]), axis = 0)
print(len(tuple_train), len(tuple_val), len(tuple_test))

print("Phasing")
#Convert tuples into Phased-tuples
phased_id2seq = {}
phased_train = []
for tup in tuple_train:
    seq = id2seq[tup[0]]
    phased_seqs = get_overlapped_chunks(seq, window, overlap)
    for phased_seq in phased_seqs:
        if len(phased_seq) >= 50 or len(phased_seqs) == 1:
            phased_id = tup[0]+'_'+str(phased_seqs.index(phased_seq))
            phased_id2seq.update({phased_id:phased_seq})
            phased_train.append(([phased_id,tup[1]]))
phased_tuple_train = np.array(phased_train)
###
phased_val = []
for tup in tuple_val:
    seq = id2seq[tup[0]]
    phased_seqs = get_overlapped_chunks(seq, window, overlap)
    for phased_seq in phased_seqs:
        if len(phased_seq) >= 50 or len(phased_seqs) == 1:
            phased_id = tup[0]+'_'+str(phased_seqs.index(phased_seq))
            phased_id2seq.update({phased_id:phased_seq})
            phased_val.append(([phased_id,tup[1]]))
phased_tuple_val = np.array(phased_val)
###
phased_test = []
phased_test_short = []
for tup in tuple_test:
    seq = id2seq[tup[0]]
    phased_seqs = get_overlapped_chunks(seq, window, overlap)
    for phased_seq in phased_seqs:
        if len(phased_seq) >= 50 or len(phased_seqs) == 1:
            phased_id = tup[0]+'_'+str(phased_seqs.index(phased_seq))
            phased_id2seq.update({phased_id:phased_seq})
            phased_test.append(([phased_id,tup[1]]))
            if len(phased_seq) <=50 and len(phased_seqs) ==1: # To get a small set of short seqs to test against
                phased_test_short.append(([phased_id, tup[1]]))
phased_tuple_test = np.array(phased_test)
phased_tuple_test_short = np.array(phased_test_short)

print("Embedding")
prot2embed = {}
for idx in phased_id2seq:
    prot2embed[idx] = to_onehot(phased_id2seq[idx], aaindex) # was juust idx before
    #print(prot2embed)
    #sys.exit("s")
    #prot2embed[idx] = to_onehot(id2seq[idx], aaindex) # change to every splitted seq
print(len(prot2embed))


#config = tf.ConfigProto(device_count = {'GPU':1})
config = tf.ConfigProto(device_count = {'CPU':1})
#config =  tf.compat.v1.ConfigProto()
#config.gpu_options.allow_growth = True
#config.gpu_options.per_process_gpu_memory_fraction = 0.3
sess = tf.Session(config=config)
session = tf.compat.v1.keras.backend.get_session()
seq, flat = get_seq_model(params)
output = Dense(1, activation='sigmoid')(flat)
model = Model(inputs=[seq], outputs=output)
model.summary()
#model = multi_gpu_model(model, 1)
model.compile(
loss='binary_crossentropy',
optimizer=Adam(),
metrics=['accuracy']) # Can stop here for indivudial fasta testing
checkpoint = ModelCheckpoint(model_file, monitor='val_accuracy', verbose = 1, save_best_only=True, mode='max')

generator = seq_Generator(phased_tuple_train[:,0], phased_tuple_train[:,1], batch_size, prot2embed, MAXLEN=MAXLEN)
test_generator = seq_Generator(phased_tuple_test[:,0], phased_tuple_test[:,1], batch_size, prot2embed, MAXLEN=MAXLEN)
short_test_generator = seq_Generator(phased_tuple_test_short[:,0], phased_tuple_test_short[:,1], batch_size, prot2embed, MAXLEN=MAXLEN)
val_generator = seq_Generator(phased_tuple_val[:,0], phased_tuple_val[:,1], batch_size, prot2embed, MAXLEN=MAXLEN)



history = model.fit_generator(generator=generator,
                  validation_data=val_generator,
                  epochs=25,
                  steps_per_epoch = 200,#tuple_train.shape[0]//batch_size,#100,
                  verbose=1,
                  callbacks=[checkpoint])

model.save('final_model_Taxa')

print("Accuracy Scoring:")
from sklearn.metrics import accuracy_score
predict = model.predict_generator(test_generator)
test_label = phased_tuple_test[:,1]
test_label = test_label.tolist()
test_label = ([int(x) for x in test_label])

predict_label = []
for score in predict:
    score = float(score[0])
    if score >= 0.5:
        predict_label.append(1)
    else:
        predict_label.append(0)
print("Accuracy_HERE")
print(len(predict))
print(accuracy_score(test_label,predict_label))
print(predict)

Correct = 0

incorrect_seqs_phased = open('incorrect_taxa_phased.fa', 'w')
incorrect_seqs = open('incorrect_taxa.fa', 'w')

incorr = []

for idx, p in enumerate(predict_label):
    if p == 0 and test_label[idx] == 0:
        Correct +=1
    elif p == 1 and test_label[idx] == 1:
        Correct +=1
    else:
        incorrect_id = phased_tuple_test[idx]
        incorrect_seqs_phased.write('>_Incorrect_' + incorrect_id[0] + '\n')
        incorrect_seqs_phased.write(phased_id2seq[incorrect_id[0]] + '\n')


        ic = incorrect_id[0]
        ic  = ic.rpartition('_')
        if ic[0] not in incorr:
            incorrect_seqs.write('>_Incorrect_' + incorrect_id[0] + '\n')
            incorr.append(ic[0])
            incorrect_seqs.write(id2seq[ic[0]] + '\n')

print("Number of Correct out")
print(Correct)



all_test = open('all_phased_test.fa', 'w')
for pt in phased_test:
    all_test.write('>'+pt[0]+'\n')
    all_test.write(phased_id2seq[pt[0]]+'\n')

all_train_phased = open('all_phased_train.fa', 'w')
for pt in phased_train:
    all_train_phased.write('>'+pt[0]+'\n')
    all_train_phased.write(phased_id2seq[pt[0]]+'\n')

all_train = open('all_train.fa', 'w')
for tup in tuple_train:
    seq = id2seq[tup[0]]
    all_train.write('>'+tup[0]+'\n')
    all_train.write(seq+'\n')

all_seqs = open('all_phased_seqs.fa', 'w')
for key,value in phased_id2seq.items():
    all_seqs.write('>'+key+'\n')
    all_seqs.write(value+'\n')

print("Short Seqs: "+str(len(phased_tuple_test_short)))
predict = model.predict_generator(short_test_generator)
test_label = phased_tuple_test_short[:,1]
test_label = test_label.tolist()
test_label = ([int(x) for x in test_label])

predict_label = []
for score in predict:
    score = float(score[0])
    if score >= 0.5:
        predict_label.append(1)
    else:
        predict_label.append(0)
print("Short Accuracy_HERE")
print(len(predict))
print(accuracy_score(test_label,predict_label))
#print(predict)


