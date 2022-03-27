import numpy as np
import pandas as pd
import os
pd.set_option('display.max_rows',None)#
pd.set_option('display.max_columns',None)#
pd.set_option('display.width',1000)#
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=1000) #

def load_data(file):
    lista = []
    records = list(open(file, "r"))
    records = records[0:]
    #print(records)
    X = []
    Y = []
    for seq in records:

        elements = seq.strip().split("\t")

        #print("elements",elements)
        X.append(elements[2:])
        #print(X)
        level = elements[0].split("\n")
        #print(level)
        classe = level[0]
        ##Y.append(ciao.index(classe))
        Y.append(classe)
    #print(X)

    X = np.array(X, dtype=float)
    #print("****************")
    Y = np.array(Y, dtype=int)
    #print("y", Y)

    return X, Y,  len(X[0])

AA = 'ACDEFGHIKLMNPQRSTVWY'
pssm_ns=[]
header = []
for i in AA:
    header.append('pssm_'+str(i) )
#print(header)
pssm_ns.append(header)

path="D:\lunwen\simplifyPSSM//"  #

path_list=os.listdir(path)
path_list.sort(key=lambda x: int(x.split('.')[0]))#

for filename in path_list:
    print('filename:', os.path.join(path, filename))
    X, labels, input_length = load_data( os.path.join(path, filename))
    #print(X)
    print(X.sum(axis=0))#
    pssm_n = []
    print(len(X))
    print("************")
    length=len(X)
    #length = 50
    for i in range(20):
        a = 0
        for j in range(int(length)):     #len(X)
            # print(X[j][i])
            a += X[j][i]
        aa=a/length
        pssm_n.append(aa)
        #print(a)
    #print(pssm_n)
    pssm_ns.append(pssm_n)
#print(pssm_ns)

def savetsv(encodings, file = 'aac_pssm.tsv'):
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

outFile = 'aac_pssm.tsv'
savetsv(pssm_ns, outFile)

