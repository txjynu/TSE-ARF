import sys, os
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import re
from collections import Counter
import numpy as np
import pandas as pd
import sys
#f = open('AAC_PSSM_ZSCALE.log', 'a')
#sys.stdout = f
#sys.stderr = f # redirect std err, if necessary

def load_data(file):
    lista = []
    records = list(open(file, "r"))
    records = records[1:]
    for seq in records:
        elements = seq.split("\t")
        level = elements[0].split("\n")
        classe = level[0]
        lista.append(classe)

    lista = set(lista)
    classes = list(lista)
    #print("class",classes)
    X = []
    Y = []
    for seq in records:
        elements = seq.strip().split("\t")
        # print("elements",elements)
        X.append(elements[1:])
        level = elements[0].split("\n")
        classe = level[0]
        ##Y.append(ciao.index(classe))
        Y.append(classe)
        #print(Y)

    X = np.array(X, dtype=float)
    #print("****************")
    Y = np.array(Y)
    #print("y", Y)

    return X, Y, len(classes), len(X[0])


zscale1 = np.array([ [0.24,  -2.32,  0.60, -0.14,  1.30], # A
		[0.84,  -1.67,  3.71,  0.18, -2.65], # C
		[3.98,   0.93,  1.93, -2.46,  0.75], # D
		[3.11,   0.26, -0.11, -0.34, -0.25], # E
		[-4.22,  1.94,  1.06,  0.54, -0.62], # F
		[2.05,  -4.06,  0.36, -0.82, -0.38], # G
		[2.47,   1.95,  0.26,  3.90,  0.09], # H
		[-3.89, -1.73, -1.71, -0.84,  0.26], # I
		[2.29,   0.89, -2.49,  1.49,  0.31], # K
		[-4.28, -1.30, -1.49, -0.72,  0.84], # L
		[-2.85, -0.22,  0.47,  1.94, -0.98], # M
		[3.05,   1.62,  1.04, -1.15,  1.61], # N
		[-1.66,  0.27,  1.84,  0.70,  2.00], # P
		[1.75,   0.50, -1.44, -1.34,  0.66], # Q
		[3.52,   2.50, -3.50,  1.99, -0.17], # R
		[2.39,  -1.07,  1.15, -1.39,  0.67], # S
		[0.75,  -2.18, -1.12, -1.46, -0.40], # T
		[-2.59, -2.64, -1.54, -0.85, -0.02], # V
		[-4.36,  3.94,  0.59,  3.44, -1.59], # W
		[-2.54,  2.44,  0.43,  0.04, -1.47] # Y
                ])

X, labels,nb_classes, input_length = load_data(  'D:\lunwen\AAC-PSSM_50.tsv')


AA = 'ACDEFGHIKLMNPQRSTVWY'
codings = []
header = []
for i in AA:
    for z in ('1', '2', '3', '4', '5'):
        header.append('acc_cpssm_'+str(i) + '.Z' + z) #
#print(header)
codings.append(header)

for i in X:
    #print(i)
    coding = []
    code=i
    #print("***************")
    ii = 0
    for i in code:
        #print(i)
        encod=[]
        for j in zscale1[ii]:
            #print(zscale1[ii])
            encod.append(i*j)
        ii=ii+1
        #print("encod",encod)
        coding.extend(encod)

    #print(coding)
    codings.append(coding)
print(codings)

#print( np.matrix(codings))
import sys

def savetsv(encodings, file = 'ZAAC-PSSM.tsv'):
	with open(file, 'w') as f:
		if encodings == 0:
			f.write('Descriptor calculation failed.')
		else:
			for i in range(len(encodings[0]) - 1):
				f.write(encodings[0][i] + '\t')
			f.write(encodings[0][-1] + '\n')
			for i in encodings[1:]:
				#f.write(i[0] + '\t')
				for j in range(0, len(i) - 1):
					f.write(str(float(i[j])) + '\t')
				f.write(str(float(i[len(i)-1])) + '\n')
	return None
outFile = 'ZAAC-CPSSM.tsv'
savetsv(codings, outFile)

