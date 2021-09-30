import numpy as np
import collections

results = np.load('../model/classifying.npy')
#phased_tuple_test = np.load('../model/phased_tuple_test.npy')

correctness = collections.defaultdict(int)

input = open('../TEST.csv','r')

data = []
for line in input:
    data.append(line)

for idx, classf in enumerate(data):
    classf = int(classf[1])
    score = results[idx]
    score = score.item()
    score = round(score)
    if score == classf:
        correctness["Correct"] +=1
    else:
        correctness["Incorrect"] +=1

print(correctness["Correct"])
print(correctness["Incorrect"])