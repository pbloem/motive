# -*- coding: utf-8 -*-

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as p
import numpy as n
import pylab
import scipy.stats as stats
import networkx as nwx
import glob
import builtins 
from matplotlib.pyplot import margins
import os.path
import json

from sklearn import svm, cross_validation, datasets


def classify(data, cls):
    y = data[0, :]
    X = data[1:, :]
    
    if(cls == 'svm'):
        model = svm.SVC(kernel='linear');
        X_train, X_test, y_train, y_test = cross_validation.train_test_split(X, y, test_size=0.2)
        
    
    # split the data 80/20

for path in sorted(glob.glob('motive.*.csv')):
    data = n.genfromtxt(path, delimiter=',')
    svmresult = classify(data, 'svm')

for path in sorted(glob.glob('motive.*.csv')):
    data = n.genfromtxt(path, delimiter=',')
    svmresult = classify(data, 'svm')

 


ax = p.subplot(111)
ax.plot(data[:, :], alpha=0.5)

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(True)
ax.spines["left"].set_visible(True)

ax.get_xaxis().set_tick_params(which='both', top='off')

ax.set_xlabel('iterations')
ax.set_ylabel('perturbation')

p.savefig('am.perturbation.png')