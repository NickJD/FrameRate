import numpy as np
import random
from utils import *
from random import seed
from random import randint
import sys
from random import shuffle
from sklearn.metrics import accuracy_score

subsample_factor = float(sys.argv[1])
params = {}
params["min_kernel"] = int(sys.argv[2])
params["n_filters"] = int(sys.argv[3])
params["pool_size"] = int(sys.argv[4])
params['dense_units'] = int(sys.argv[5])
print("Params: "+str(params))

model_file = 'saved_models/model_s' + str(subsample_factor) + "_k" + str(params["min_kernel"]) + "_f" + str(params["n_filters"]) + \
                  "_p" + str(params["pool_size"]) + "_d" + str(params["dense_units"])
print(model_file)

MAXLEN = 75
batch_size = 2000
window = 75
overlap = 25
id2seq = {}
aaletters = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
aaindex = {'K': 0, 'H': 1, 'T': 2, 'V': 3, 'Q': 4, 'L': 5, 'D': 6, 'R': 7, 'A': 8, 'G': 9, 'Y': 10, 'F': 11,
           'S': 12, 'N': 13, 'P': 14, 'C': 15, 'M': 16, 'I': 17, 'W': 18, 'E': 19}

family2count = {}
geneid2seq = {}
with open("../eggnog/ensembl_cds_genes_combined_cdhit_eggnogmapped_COGed_6AAs.fasta", 'r') as f:
    for line in f:
        if ">" == line[0]:
            idx, label = line.strip().split(",")
            cog, gene, frame = idx.split("_")
            geneid = cog+"_"+gene
            seq = next(f).strip()
            
            if label == "1":
                if cog not in family2count:
                    family2count[cog] = 0
                family2count[cog] += 1
            
            
#             if "*" in seq:
#                 if label == "1" and seq[-1] != "*": 
#                     print("* in positive seq ", geneid)
                
print(len(family2count))

# 5000 each
thres = [0, 10, 100, 1000, 5000, 20000]
group2families = {}
print("l, h, No.Families")
for i in range(1, len(thres)):
    group2families[thres[i]] = [k for k, v in family2count.items() if v > thres[i-1] and v <= thres[i]]
    print(thres[i-1], thres[i], len(group2families[thres[i]]))
    
def split_cogs(existing):
    families = set()
    print("No.GenesinGroup, No.Families")
    for i in range(1, len(thres)):
        total = 0
        while total < 5000: 
            sample = random.choice(group2families[thres[i]])
            if existing == None or (existing != None and sample not in existing):
                families.add(sample)
                total += family2count[sample]
        print(total, len(families))
    return families

validFamilies = split_cogs(None)
testFamilies = split_cogs(validFamilies)
trainFamilies = family2count.keys() - validFamilies - testFamilies

if len(testFamilies & validFamilies & trainFamilies) == 0:
    print("there is no overlap")
    
def get_overlapped_chunks(textin, chunksize, overlapsize):  # Maybe remove last chunk as it will be covered and will be short???
    phases = []
    # if len(textin) > chunksize:
    end_flag = False
    for phase in [ textin[a:a+chunksize] for a in range(0,len(textin), chunksize-overlapsize) ]: # 75 - 25 = 50
        if len(phase) < chunksize and end_flag != True:
            end_flag = True
            phases.append(textin[-chunksize:])
        elif len(phase) == chunksize:
            phases.append(phase)
    return phases

x = {"train_pos":[], "train_neg":[], "test_pos":[],"test_neg":[], "val_pos":[], "val_neg":[]}
y = {"train_pos":[], "train_neg":[], "test_pos":[],"test_neg":[], "val_pos":[], "val_neg":[]}

window = 75
overlap = 25

count_negWithStop = 0
count_train = 0
with open("../eggnog/ensembl_cds_genes_combined_cdhit_eggnogmapped_COGed_6AAs.fasta", 'r') as f:
    for line in f:
        if ">" == line[0]:
            idx, label = line.strip().split(",")
            cog, gene, frame = idx.split("_")
            geneid = cog+"_"+gene
            seq = next(f).strip()
            
            if label == "1" and seq[-1] != "*":
                seq = seq[:-1]
            if "*" not in seq: # filter negative frames with "*"
                
                flag = "pos" if label == "1" else "neg"
                if cog in trainFamilies:
                    count_train += 1
                    if count_train > subsample_factor:
                        continue
                    group = "train"
                elif cog in validFamilies:
                    group = "val"
                else:
                    group = "test"
                key = group + "_" + flag
                
                chunks = get_overlapped_chunks(seq, window, overlap)
                for chunk in chunks:
                    onehot = to_onehot(chunk, aaindex)
                    x[key].append(onehot)
            else:
                count_negWithStop += 1            
            
#             if "*" in seq:
#                 if label == "1" and seq[-1] != "*": 
#                     print("* in positive seq ", geneid)
for key in x:
    x[key] = np.array(x[key])
    np.random.shuffle(x[key])
    print(key, x[key].shape)
    
for split in ["train", "test", "val"]:
    x[split+"_neg"] = x[split+"_neg"][:len(x[split+"_pos"]), :]    

for key in x:
    if "pos" in key:
        y[key] = np.ones(len(x[key]))
    else:
        y[key] = np.zeros(len(x[key]))
    print(key, x[key].shape, y[key].shape)
    X, Y = {}, {} 
    
for split in ["train", "test", "val"]:
    X[split] = np.concatenate((x[split+"_pos"], x[split+"_neg"]), axis = 0)
    print(X[split].shape)
    Y[split] = np.concatenate((y[split+"_pos"], y[split+"_neg"]), axis = 0)
    print(Y[split].shape)
    
runs = 1 
accuracies = []
for i in range(runs):
    seq, flat = get_seq_model(params)
    output = Dense(1, activation='sigmoid')(flat)
    model = Model(inputs=[seq], outputs=output)
    model.summary()
    model.compile(
    loss='binary_crossentropy',
    optimizer=Adam(),
    metrics=['accuracy']) # Can stop here for indivudial fasta testing
    
    checkpoint = ModelCheckpoint(model_file, save_weights_only=True, monitor='val_accuracy', verbose = 1, save_best_only=True, mode='max')

    history = model.fit(X["train"], Y["train"], validation_data = (X["val"], Y["val"]), validation_steps = X["val"].shape[0]//batch_size,
                      epochs=5,
                      steps_per_epoch = X["train"].shape[0]//batch_size,
                      verbose=1,
                      callbacks=[checkpoint], workers = 20,
                      use_multiprocessing=True,
                      max_queue_size=100)
    
    model.load_weights(model_file)
    print("weights loaded")
    model.save(model_file+".h5")
    print("full model saved")
    del model 
    print("model deleted")
    model = load_model(model_file+".h5")
    print("model loaded")
    predict = model.predict(X["test"])
    # test_label = ([int(x) for x in test_label])

    predict_label = []
    for score in predict:
        score = float(score[0])
        if score >= 0.5:
            predict_label.append(1)
        else:
            predict_label.append(0)

    print(len(predict))
    acc = accuracy_score(Y["test"], predict_label)
    print("Accuracy", acc)
    accuracies.append(acc)

print("Accuracy from each run is ", accuracies)
print("Mean is ", np.mean(accuracies))
    
