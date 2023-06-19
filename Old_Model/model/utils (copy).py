from keras.callbacks import ModelCheckpoint, EarlyStopping, CSVLogger
import numpy as np
from keras.models import Model, load_model
from keras.layers import (
    Input, Dense, Embedding, Conv1D, Flatten, Concatenate,
    MaxPooling1D, Dropout, Dot, LeakyReLU
)
from keras import backend as K
from tensorflow.keras.optimizers import Adam, RMSprop
from sklearn.metrics import roc_auc_score, average_precision_score
from keras.utils import np_utils
from tensorflow.keras.utils import Sequence #, multi_gpu_model
import scipy.stats as ss

def repeat_to_length(s, wanted):
    return (s * (wanted//len(s) + 1))[:wanted]

def to_onehot(seq, aaindex, MAXLEN = 75, repeat=False):
    onehot = np.zeros((MAXLEN, 20), dtype=np.int32) # after change aindex I can change to 20 here
    #original_len = min(MAXLEN, len(seq))
    if repeat == True:
        seq = repeat_to_length(seq, MAXLEN)
    for i in range(len(seq)): # used to be original_len
        onehot[i, aaindex.get(seq[i])] = 1
    #onehot[len(seq):, 0] = 1 # used to be original_len
    return onehot
    
def get_seq_model(params,MAXLEN = 75):
    seq = Input(shape=(MAXLEN, 20), dtype=np.float32)
    kernels = range(params['min_kernel'], 48, 4)
    nets = []
    for i in range(len(kernels)):
        conv = Conv1D(
            filters=params['n_filters'],
            kernel_size=kernels[i],
            padding='valid',
            kernel_initializer='glorot_normal')(seq)
        conv_dropout = Dropout(0.5)(conv)
        pool = MaxPooling1D(pool_size=params['pool_size'])(conv_dropout)
        
        conv2 = Conv1D(
            filters=params['n_filters'],
            kernel_size=kernels[i]//2,
            padding='valid',
            kernel_initializer='glorot_normal')(pool)
        conv_dropout2 = Dropout(0.5)(conv2)
        pool2 = MaxPooling1D(pool_size=params['pool_size'])(conv_dropout2)
        flat = Flatten()(pool2)
        
        nets.append(flat)

    net = Concatenate(axis=1)(nets)
    dense_seq = Dense(params['dense_units'])(net)
    activ_seq = LeakyReLU(alpha=0.1)(dense_seq)
    dropout_seq = Dropout(0.5)(activ_seq)

    return seq, dropout_seq
