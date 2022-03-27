from sklearn.model_selection import KFold
from sklearn.metrics import roc_curve, auc, accuracy_score, precision_score, recall_score, f1_score,matthews_corrcoef
from sklearn.metrics import confusion_matrix as cm
import matplotlib.pyplot as plt
import pandas as pd
import time
import random
import numpy
from sklearn.ensemble import RandomForestClassifier
from numpy import *
from sklearn import *
from pandas import *
from scipy import stats

def SFLA_RF(x_traincv, y_traincv, x_testcv, y_testcv,n_estimators=False, max_features=False,max_depth=False,min_samples_split=False):
 print(min_samples_split)
 rf = RandomForestClassifier(n_estimators=int(n_estimators),max_features=int(max_features),max_depth=int(max_depth),min_samples_split=int(min_samples_split), random_state=920).fit(x_traincv,y_traincv)
 y_score = rf.predict_proba(x_testcv)[:, 1]
 fpr, tpr, threshold = metrics.roc_curve(y_testcv, y_score, pos_label=1)
 roc_auc = metrics.auc(fpr, tpr)

 return roc_auc


def SFLA_RF_CV(x_train, y_train, n, n_estimators=False, max_features=False, max_depth=False,min_samples_split=False):
 '''
 n: number of splits for k-fold
 '''
 KF = KFold(n_splits=n, shuffle=True, random_state=920)
 f = []
 for train_indexcv, test_indexcv in KF.split(x_train):
  x_traincv, x_testcv = x_train.iloc[train_indexcv][:], x_train.iloc[test_indexcv][:]
  y_traincv, y_testcv = y_train.iloc[train_indexcv][:], y_train.iloc[test_indexcv][:]
  fq = SFLA_RF(x_traincv, y_traincv, x_testcv, y_testcv,n_estimators, max_features, max_depth,min_samples_split)
  f.append(fq)
 f = mean(f)
 return f


def SFLA_RF_rf(num_parameter, num_global, num_local, m, n, q, n1,rangeN_estimators, rangeMax_features, rangeMax_depth,rangeMin_samples_split, x_train, y_train):
 '''
  num_parameter: int, number of parameter to optimize
 num_global: int, the maximum number of global iterations
 num_local: int, the maximum number of local iterations
 m : int, the number of memeplexes
 n : int, the number of frogs in each memeplex
 q : int, the number of frogs in submemeplex
 n1:  number of splits for cross validation for inner loop
 rangeC: list, float, range of parameter C,eg.[10**-2, 10**2]
 rangeGamma: list, float, range of parameter Gamma,eg.[10**-6, 1]
 x_train: feature
 y_train: lable
 '''

 # --- Step 0--Initialize parameters ---#
 sizeN_estimators = 2
 sizeMax_features = 2
 sizeMax_depth = 2
 sizeMin_samples_split = 2
 max_step = [(rangeN_estimators[1] - rangeN_estimators[0]) / sizeN_estimators, (rangeMax_features[1] - rangeMax_features[0]) / sizeMax_features,
             (rangeMax_depth[1] - rangeMax_depth[0]) / sizeMax_depth,
             (rangeMin_samples_split[1] - rangeMin_samples_split[0]) / sizeMin_samples_split]  # maximum step size

 # --- Step 1--Generate initial population ---#
 frogN_estimators = 10 ** random.uniform(log10(rangeN_estimators[0]), log10(rangeN_estimators[1]), m * n)
 frogMax_features = 10 ** random.uniform(log10(rangeMax_features[0]), log10(rangeMax_features[1]), m * n)
 frogMax_depth = random.randint(rangeMax_depth[0], rangeMax_depth[1], m * n)
 frogMin_samples_split = random.uniform(rangeMin_samples_split[0], rangeMin_samples_split[1], m * n)
 frog = c_[frogN_estimators, frogMax_features, frogMax_depth, frogMin_samples_split]

 # Compute the performance value for each frog on validation data #
 KF = KFold(n_splits=n1, shuffle=True, random_state=920)
 f = zeros((m * n, n1))
 j = 0
 for train_indexcv, test_indexcv in KF.split(x_train):
  x_traincv, x_testcv = x_train.iloc[train_indexcv][:], x_train.iloc[test_indexcv][:]
  y_traincv, y_testcv = y_train.iloc[train_indexcv][:], y_train.iloc[test_indexcv][:]
  for i in range(m * n):
   f[i, j] = SFLA_RF(x_traincv, y_traincv, x_testcv, y_testcv, frog[i,0],frog[i,1],frog[i,2],frog[i,3])
  j += 1
 f = f.mean(axis=1)
 f_parameter = c_[f, frog]

 # --- Step 2--Rank frogs ---#
 f_parameter = f_parameter[argsort(f_parameter[:, 0])[::-1]]

 #######--- Global search start---######
 i_global = 0
 flag = 0
 fBest_iteration = f_parameter[0, 0]
 weights = [2 * (n + 1 - j) / (n * (n + 1)) for j in range(1, n + 1)]  # weights of ranked frogs in each memeplex
 while i_global < num_global:
  frog_gb = f_parameter[0, 0]  # mark the global best frog
  # --- Step 3--Partition frogs into memeplexes ---#
  memeplexes = zeros((m, n, num_parameter + 1))  # [memeplexes, frog in memeplex,[f,C,Gamma] ]
  for i in range(m):
   memeplexes[i] = f_parameter[linspace(i, m * n + i, num=n, endpoint=False, dtype=int)]

  #######--- Local search start---######
  # --- Step 4--Memetic evolution within each memeplex ---#
  im = 0  # the number of memeplexes that have been optimized
  while im < m:
   i_local = 0  # counts the number of local evolutionary steps in each memeplex
   while i_local < num_local:

    # --- Construct a submemeplex ---#
    rValue = random.random(n) * weights  # random value with probability weights
    subindex = sort(argsort(rValue)[::-1][0:q])  # index of selected frogs in memeplex
    submemeplex = memeplexes[im][subindex]  # form submemeplex

    # --- Improve the worst frog's position ---#
    # Learn from local best Pb #
    Pb = submemeplex[0]  # mark the best frog in submemeplex
    Pw = submemeplex[q - 1]  # mark the worst frog in memeplex
    S = (Pb - Pw)[1:] * (Pb - Pw)[0]
    Uq = Pw[1:] + S
    # Check feasible space and the performance #
    if (rangeN_estimators[0] <= Uq[0] <=rangeN_estimators[1]) and (rangeMax_features[0] <= Uq[1] <=rangeMax_features[1])\
            and(rangeMax_depth[0] <= Uq[2] <= rangeMax_depth[1]) and (rangeMin_samples_split[0] <= Uq[3] <= rangeMin_samples_split[1]):  # check feasible space
     fq = SFLA_RF_CV(x_train, y_train, n1, Uq[0], Uq[1],Uq[2],Uq[3])
     if fq < Pw[0]:  # if no improvement of performance,learn from global best randomly #
      S = random.random(num_parameter) * (frog_gb - Pw)[1:]
      for i in range(num_parameter):
       if S[i] > 0:
        S[i] = min(S[i], max_step[i])
       else:
        S[i] = min(S[i], -max_step[i])
      Uq = Pw[1:] + S
      if (rangeN_estimators[0] <= Uq[0] <= rangeN_estimators[1]) and (
              rangeMax_features[0] <= Uq[1] <= rangeMax_features[1]) \
              and (rangeMax_depth[0] <= Uq[2] <= rangeMax_depth[1]) and (
              rangeMin_samples_split[0] <= Uq[3] <= rangeMin_samples_split[1]):  # check feasible space
       fq = SFLA_RF_CV(x_train, y_train, n1, Uq[0],Uq[1],Uq[2],Uq[3])
       if fq < Pw[0]:  # if no improvement of performance, randomly generate a new frog
        Uq = [10**random.uniform(log10(rangeN_estimators[0]),log10(rangeN_estimators[1])),10**random.uniform(log10(rangeMax_features[0]),log10(rangeMax_features[1])),\
              10 ** random.uniform(log10(rangeMax_depth[0]), log10(rangeMax_depth[1])),
              10 ** random.uniform(log10(rangeMin_samples_split[0]), log10(rangeMin_samples_split[1]))]
        fq = SFLA_RF_CV(x_train, y_train, n1, Uq[0],Uq[1],Uq[2],Uq[3])
      else:  # if not in the feasible space, randomly generate a new frog
       Uq = [10**random.uniform(log10(rangeN_estimators[0]),log10(rangeN_estimators[1])),10**random.uniform(log10(rangeMax_features[0]),log10(rangeMax_features[1])),\
             10 ** random.uniform(log10(rangeMax_depth[0]), log10(rangeMax_depth[1])),
             10 ** random.uniform(log10(rangeMin_samples_split[0]), log10(rangeMin_samples_split[1]))]
       fq = SFLA_RF_CV(x_train, y_train, n1, Uq[0], Uq[1],Uq[2],Uq[3])
    else:  # if not in the feasible space, learn from global best randomly
     S = random.random(num_parameter) * (frog_gb - Pw)[1:]
     for i in range(num_parameter):
      if S[i] > 0:
       S[i] = min(S[i], max_step[i])
      else:
       S[i] = min(S[i], -max_step[i])
     Uq = Pw[1:] + S
     if  (rangeN_estimators[0] <= Uq[0] <=rangeN_estimators[1]) and (rangeMax_features[0] <= Uq[1] <=rangeMax_features[1])\
            and(rangeMax_depth[0] <= Uq[2] <= rangeMax_depth[1]) and (rangeMin_samples_split[0] <= Uq[3] <= rangeMin_samples_split[1]):  # check feasible space
      fq = SFLA_RF_CV(x_train, y_train, n1, Uq[0], Uq[1],Uq[2],Uq[3])
      if fq < Pw[0]:  # if no improvement of performance, randomly generate a new frog
       Uq = [10**random.uniform(log10(rangeN_estimators[0]),log10(rangeN_estimators[1])),10**random.uniform(log10(rangeMax_features[0]),log10(rangeMax_features[1])),
             10 ** random.uniform(log10(rangeMax_depth[0]), log10(rangeMax_depth[1])),
             10 ** random.uniform(log10(rangeMin_samples_split[0]), log10(rangeMin_samples_split[1]))]
       fq = SFLA_RF_CV(x_train, y_train, n1, Uq[0],  Uq[1],Uq[2],Uq[3])
     else:  # if not in the feasible space, randomly generate a new frog
      Uq = [10**random.uniform(log10(rangeN_estimators[0]),log10(rangeN_estimators[1])),10**random.uniform(log10(rangeMax_features[0]),log10(rangeMax_features[1])),
            10 ** random.uniform(log10(rangeMax_depth[0]), log10(rangeMax_depth[1])),
            10 ** random.uniform(log10(rangeMin_samples_split[0]), log10(rangeMin_samples_split[1]))]
      fq = SFLA_RF_CV(x_train, y_train, n1, Uq[0],  Uq[1],Uq[2],Uq[3])

    # --- Upgrade the memeplex ---#
    memeplexes[im][subindex[q - 1]] = r_[fq, Uq]
    memeplexes[im] = memeplexes[im][argsort(memeplexes[im][:, 0])[::-1]]

    i_local += 1

   im += 1
  #######--- Local search end---######

  # --- Step 5--Shuffle memeplexes ---#
  f_parameter = memeplexes.reshape(m * n, num_parameter + 1)
  f_parameter = f_parameter[argsort(f_parameter[:, 0])[::-1]]

  i_global += 1

  # --- Step 6--Check convergence ---#
  if f_parameter[0, 0] > 0.999:
   print('The program was terminated because it reached the optimization goal with f = %.3f' % f_parameter[0, 0])
   break
  if abs(frog_gb - f_parameter[0, 0]) < 10 ** -4:
   flag += 1
  if flag > 30:
   break
  fBest_iteration = r_[fBest_iteration, f_parameter[0, 0]]

  #######--- Global search end---######

 return (f_parameter[0], fBest_iteration)

# --- Ensemble ---#
def OptimizeRF_SFLA_CV(x, y, n_splits, num_parameter, num_global, num_local, m, n, q, n1,rangeN_estimators, rangeMax_features, rangeMax_depth,rangeMin_samples_split):
 ##---Classification with n-fold cross-validation---##
 # --- x is feature, y is lable, n is number of fold
 # ---  define K-fold cross validation ---#
 KF = KFold(n_splits, shuffle=True, random_state=920)
 y_score = []
 y_test = []
 for train_index, test_index in KF.split(x):
  # ---  Seperate traing set and test set ---#
  x_train, x_test = x.iloc[train_index][:], x.iloc[test_index][:]
  y_train = y.iloc[train_index][:]

  # ---  Fill NaN age ---#
  #x_train[isnan(x_train)] = 0
  #x_test[isnan(x_test)] = 0
  #print(y_train.shape)

  ##---  optimize RF with SFLA---##
  x_train = pd.DataFrame(x_train)
  y_train = pd.Series(y_train)
  f_parameter, fBest_iteration = SFLA_RF_rf(num_parameter, num_global, num_local, m, n, q, n1,rangeN_estimators, rangeMax_features,rangeMax_depth,rangeMin_samples_split
                                         , x_train, y_train)
  # f_parameter: list, [bestAUC,bestC,bestGamma,bestDegree,bestCoef0]   fBest_iteration: bestAUC in each iteration

  ##---  creat and train the model ---##
  rf =  RandomForestClassifier(n_estimators=int(f_parameter[1]),max_features=int(f_parameter[2]),max_depth=int(f_parameter[3]),min_samples_split=int(f_parameter[4]), random_state=920)

  print(f_parameter)
  rf.fit(x_train, y_train)
  # Plot ROC and calculate AUC
  y_score.extend([x[1] for x in rf.predict_proba(x_test).tolist()])
  y_test.extend(y[test_index].tolist())
  fpr, tpr, threshold = roc_curve(y[test_index], rf.predict_proba(x_test)[:, 1], pos_label=1)
  roc_auc = auc(fpr, tpr)
  print('AUC:', roc_auc)

import time

# import datas/et1
x = pd.read_csv('D:\lunwen\T3SE.csv',index_col=0)
y=pd.read_csv('D:\lunwen\T3SS.csv',index_col=0)
y=y.iloc[:,0]
print(x.shape)
print(y.shape)

start = time.process_time()
n_splits = 5  #10 number of splits for outer loop
num_parameter = 4  # number of parameter to optimize
num_global = 5  # the maximum number of global iterations
num_local = 5  # the maximum number of local iterations
m = 4  # the number of memeplexes
n = 8  # the number of frogs in each memeplex
q = 5  # the number of frogs in submemeplex
n1 = 5  #10 number of splits for inner loop
rangeN_estimator = [10, 500]  # 
rangeMax_features = [2, 240]  #
rangeMax_depth = [2, 100]  #
rangeMin_samples_split = [2, 100]  # 
rf, a, p, r, f1score, roc_auc, y_pred, y_score, fpr, tpr = OptimizeRF_SFLA_CV(x, y, n_splits, num_parameter,
                                                                                   num_global, num_local, m, n, q, n1,
                                                                                   rangeN_estimator, rangeMax_features,
                                                                                   rangeMax_depth,
                                                                                   rangeMin_samples_split)
end = time.process_time()
print('OptimizeRF_SFLA_CV algorithm takes ' + str(end - start) + 'seconds.\n')

import csv

fname = 'ARF_result.csv'

# header
tmp = [['Id', 'Prediction', 'Probability']]

# add ID numbers for each Y
for (i, y) in enumerate(y_pred):
 tmp2 = [(i + 1), y_pred, y_score]
 tmp.append(tmp2)

# write CSV file
with open(fname, 'w') as f:
 writer = csv.writer(f)
 writer.writerows(tmp)


