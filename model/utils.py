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

def get_params():
    pi = 11
    params = {}
    if pi != -1:
        max_kernels = [17, 33, 65]
        nb_filters = [8, 16]
        dense_units = [8, 16, 32]
        pool_sizes = [10, 20]#200
        params['max_kernel'] = max_kernels[pi % len(max_kernels)]
        pi //= len(max_kernels)
        params['nb_filters'] = nb_filters[pi % len(nb_filters)]
        pi //= len(nb_filters)
        params['pool_size'] = pool_sizes[pi % len(pool_sizes)]
        pi //= len(pool_sizes)
        params['dense_units'] = dense_units[pi % len(dense_units)]
        pi //= len(dense_units)
    print('Params:', params)
    return params


def split_train(data, split1, split2):
    train = data[:int(split1*len(data))]
    val = data[int(split1*len(data)):int(split2*len(data))]
    test = data[int(split2*len(data)):]
    return train, val, test

def repeat_to_length(s, wanted):
    return (s * (wanted//len(s) + 1))[:wanted]


#MAXLEN used to be 22 but moved to 26 because new AAs from Swissprot
def to_onehot(seq, aaindex, MAXLEN = 75, repeat=True):
    onehot = np.zeros((MAXLEN, 21), dtype=np.int32)
    #original_len = min(MAXLEN, len(seq))
    if len(seq) < 75:
        print("E")
    if repeat == True:
        seq = repeat_to_length(seq, MAXLEN)
    for i in range(len(seq)): # used to be original_len
        onehot[i, aaindex.get(seq[i])] = 1
    #onehot[len(seq):, 0] = 1 # used to be original_len
    return onehot


class seq_Generator(Sequence):
    def __init__(self, x_set, y_set, batch_size, prot2embed, MAXLEN = 75):
        self.x, self.y = x_set, y_set
        self.batch_size = batch_size
        self.nbatch = int(np.ceil(len(self.x) / float(self.batch_size)))
        self.length = len(self.x)
        self.prot2embed = prot2embed
        self.MAXLEN = MAXLEN

    def __len__(self):
        return self.nbatch

    def __getitem__(self, idx):
        start = idx * self.batch_size
        batch_len = min(self.batch_size, (self.length)-start)
        X_batch = np.empty((batch_len, self.MAXLEN,21), dtype=np.float32)
        y_batch = np.empty(batch_len, dtype=np.float32)

        for ids in range(start, min((idx + 1) * self.batch_size, self.length)):
            array1 = self.prot2embed[self.x[ids]]
            X_batch[ids-start,:,:] = array1
            y_batch[ids-start] = self.y[ids]
        return [X_batch], y_batch
    
def get_seq_model(params,MAXLEN = 75):
    seq = Input(shape=(MAXLEN, 21), dtype=np.float32)
    kernels = range(8, params['max_kernel'], 8)
    nets = []
    for i in range(len(kernels)):
        conv = Conv1D(
            filters=params['nb_filters'],
            kernel_size=kernels[i],
            padding='valid',
            kernel_initializer='glorot_normal')(seq)
        conv_dropout = Dropout(0.5)(conv)
        pool = MaxPooling1D(pool_size=params['pool_size'])(conv_dropout)
        flat = Flatten()(pool)
        nets.append(flat)

    net = Concatenate(axis=1)(nets)
    dense_seq = Dense(params['dense_units'])(net)
    activ_seq = LeakyReLU(alpha=0.1)(dense_seq)
    dropout_seq = Dropout(0.5)(activ_seq)

    return seq, dropout_seq
