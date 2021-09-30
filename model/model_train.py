from tensorflow.python.client import device_lib
import tensorflow.compat.v1 as tf
import random
from utils import *
import sys


params = get_params()
params["pool_size"] = 10
#sys.exit(params)

MAXLEN = 75
batch_size = 100
phasing = 75
overlapped = 50
positives, negatives = [], []
id2seq = {}
aaletters = set()
aaletters = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'}
aaindex = dict()
model_file = "model.hdf5"


print("Setup:")
print(device_lib.list_local_devices())
print("Num GPUs Available: ", len(tf.config.experimental.list_physical_devices('GPU')))



def get_overlapped_chunks(textin, chunksize, overlapsize):  # Maybe remove last chunk as it will be covered and will be short???
    return [ textin[a:a+chunksize] for a in range(0,len(textin), chunksize-overlapsize)]





with open('../Extended_CoDing_Sequences_For_Training.csv', 'r') as f: # Training Data
    count = 0
    for line in f:
        count +=1
        items = line.strip().split(',')
        if count >= 1000000:
            break
        if aaletters.issuperset(items[3]):
            if len(items[3]) >= 50: # Minimum sequence length
                if items[2] == '1':
                    #if len(positives) < 100000:
                    positives.append((items[0],1))
                    id2seq[items[0]] = items[3]
                elif items[2] == '0': # randomly - using seed - pick one of 5 negatives for each positive
                    #if len(negatives) < 100000:
                    negatives.append((items[0], 0))
                    id2seq[items[0]] = items[3]


aaletters = list(aaletters)
for i in range(len(aaletters)):
    aaindex[aaletters[i]] = i + 1

print("Positive-Negative-IDs-AA")
print(len(positives), len(negatives), len(id2seq), len(aaletters))
shuffled_positives = random.shuffle(positives) # Extra shuffling
train_pos, val_pos, test_pos = split_train(positives, 0.8, 0.9)
train_neg, val_neg, test_neg = split_train(negatives, 0.8, 0.9)
print(len(train_pos), len(val_pos), len(test_pos))
print(len(train_neg), len(val_neg), len(test_neg))









#train_pos = np.repeat(np.array(train_pos), len(train_neg)//len(train_pos), axis = 0)
#tuple_train = np.concatenate((train_pos, np.array(train_neg)), axis=0)
tuple_train = np.concatenate((train_pos, train_neg[:len(train_pos)]), axis=0)
np.random.shuffle(tuple_train)
tuple_val = np.concatenate((val_pos, val_neg[:len(val_pos)]), axis = 0)
tuple_test = np.concatenate((test_pos, test_neg[:len(test_pos)]), axis = 0)
print(len(tuple_train), len(tuple_val), len(tuple_test))

#Convert tuples into Phased-tuples
phased_id2seq = {}
phased_train = []
for tup in tuple_train:
    seq = id2seq[tup[0]]
    phased_seqs = get_overlapped_chunks(seq, phasing, overlapped)
    for phased_seq in phased_seqs:
        if len(phased_seq) >= 50:
            phased_id = tup[0]+'_'+str(phased_seqs.index(phased_seq))
            phased_id2seq.update({phased_id:phased_seq})
            phased_train.append(([phased_id,tup[1]]))
phased_tuple_train = np.array(phased_train)
###
phased_val = []
for tup in tuple_val:
    seq = id2seq[tup[0]]
    phased_seqs = get_overlapped_chunks(seq, phasing, overlapped)
    for phased_seq in phased_seqs:
        if len(phased_seq) >= 50:
            phased_id = tup[0]+'_'+str(phased_seqs.index(phased_seq))
            phased_id2seq.update({phased_id:phased_seq})
            phased_val.append(([phased_id,tup[1]]))
phased_tuple_val = np.array(phased_val)
###
phased_test = []
for tup in tuple_test:
    seq = id2seq[tup[0]]
    phased_seqs = get_overlapped_chunks(seq, phasing, overlapped)
    for phased_seq in phased_seqs:
        if len(phased_seq) >= 50:
            phased_id = tup[0]+'_'+str(phased_seqs.index(phased_seq))
            phased_id2seq.update({phased_id:phased_seq})
            phased_test.append(([phased_id,tup[1]]))
phased_tuple_test = np.array(phased_test)


prot2embed = {}
for idx in phased_id2seq:
    prot2embed[idx] = to_onehot(phased_id2seq[idx], aaindex) # was juust idx before
    #print(prot2embed)
    #sys.exit("s")
    #prot2embed[idx] = to_onehot(id2seq[idx], aaindex) # change to every splitted seq
print(len(prot2embed))


config = tf.ConfigProto(device_count = {'GPU':1})
#config =  tf.compat.v1.ConfigProto()
config.gpu_options.allow_growth = True
config.gpu_options.per_process_gpu_memory_fraction = 0.3
sess = tf.Session(config=config)

session = tf.compat.v1.keras.backend.get_session()
#K.set_session(session)
#K.set_session(sess)


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


#i = phased_tuple_train[:,0]

#e = phased_tuple_train[:,1]


generator = seq_Generator(phased_tuple_train[:,0], phased_tuple_train[:,1], batch_size, prot2embed, MAXLEN=MAXLEN)
test_generator = seq_Generator(phased_tuple_test[:,0], phased_tuple_test[:,1], batch_size, prot2embed, MAXLEN=MAXLEN)
val_generator = seq_Generator(phased_tuple_val[:,0], phased_tuple_val[:,1], batch_size, prot2embed, MAXLEN=MAXLEN)

history = model.fit_generator(generator=generator,
                  validation_data=val_generator,
                  epochs=10,
                  steps_per_epoch = 10000,#tuple_train.shape[0]//batch_size,#100,
                  verbose=1,
                  callbacks=[checkpoint])

predict = model.predict_generator(test_generator)


test_label = phased_tuple_test[:,1]
test_label = test_label.tolist()
test_label = ([int(x) for x in test_label])

from sklearn.metrics import accuracy_score

predict_label = []

for score in predict:
    score = float(score[0])
    if score >= 0.5:
        predict_label.append(1)
    else:
        predict_label.append(0)
print("Accuracy_HERE")
print(accuracy_score(test_label,predict_label))
print(predict)


#print(predict)

# for idx, val in enumerate(phased_tuple_test):
#     print(val+'_'+predict[idx])


#test_generator = seq_Generator(phased_tuple_test[:,0], phased_tuple_test[:,1], batch_size, prot2embed, MAXLEN=MAXLEN)


#classifying = model.predict_generator(test_generator)



#print(classifying)
#np.save('classifying.npy',classifying)

#print(phased_tuple_test)
#np.save('phased_tuple_test.npy',phased_tuple_test)
