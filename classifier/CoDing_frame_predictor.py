import keras
import random
mymodle = keras.models.load_model("../model/model.hdf5")
from keras.callbacks import ModelCheckpoint, EarlyStopping, CSVLogger
import numpy as np
from keras.models import Model, load_model
from keras.layers import (
    Input, Dense, Embedding, Conv1D, Flatten, Concatenate,
    MaxPooling1D, Dropout, Dot, LeakyReLU
)
from keras import backend as K
#from tensorflow.compat.v1.keras import backend as K
from tensorflow.keras.optimizers import Adam, RMSprop
from sklearn.metrics import roc_auc_score, average_precision_score
from keras.utils import np_utils
from tensorflow.keras.utils import Sequence

import scipy.stats as ss
#import tensorflow as tf
import tensorflow.compat.v1 as tf
from utils import *
import sys
import re


import os
#os.environ['CUDA_VISIBLE_DEVICES'] = "1"


params = get_params()
MAXLEN = 100
batch_size = 500
#weights_file = "frame.hdf5"

phasing = 75
overlapped = 50

positives = []
negatives = []
id2seq = {}
aaletters = set()

from tensorflow.python.client import device_lib
print(device_lib.list_local_devices())

print("Num GPUs Available: ", len(tf.config.experimental.list_physical_devices('GPU')))
aaletters = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'}

#sys.exit('Stopped_Here')

pos = []
neg = []



with open('../Extended_CoDing_Sequences_For_Training_test.csv', 'r') as f:
#with open('./swiss_alt.fasta','r') as f:
    count = 0
    for line in f:
        count +=1
        items = line.strip().split(',')
        if aaletters.issuperset(items[3]):

            if len(items[3]) > 6: # dont allo short seqqs
                id2seq[items[0]] = items[3][:100]



aaindex = dict()
aaletters = list(aaletters)
for i in range(len(aaletters)):
    aaindex[aaletters[i]] = i + 1





prot2embed = {}
for idx in id2seq:
    prot2embed[idx] = to_onehot(id2seq[idx], aaindex) 
    
print("PROT")
print(prot2embed[idx])
print("END")
print(aaindex)
sys.exit("exited")

print(len(prot2embed))


tuple_run = {x: 1 for x in id2seq}
tuple_run = list(tuple_run.items())
tuple_run = np.array(tuple_run)

test_generator = seq_Generator(tuple_run[:,0], tuple_run[:,1], batch_size, prot2embed, MAXLEN=MAXLEN)


classifying = mymodle.predict_generator(test_generator)


print(tuple_run)
print("HERE")
print(classifying)

class_2 = mymodle.predict_generator(test_generator)

print(class_2)




import collections

correctness = collections.defaultdict(int)

input = open('../Extended_CoDing_Sequences_For_Training_test.csv','r')

data = {}
for line in input:
    id = line.split(',')[0]
    classif = line.split(',')[2]
    data.update({id:classif})

for idx, classf in enumerate(tuple_run):
    classf = classf[0]
    classf = int(data[classf])
    score = classifying[idx]
    score = score.item()
    score = round(score)
    if score == classf:
        correctness["Correct"] +=1
    else:
        correctness["Incorrect"] +=1

print(correctness["Correct"])
print(correctness["Incorrect"])

# confusion = tf.confusion_matrix(labels=y_, predictions=y, num_classes=num_classes)
# print(confusion)
