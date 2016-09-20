# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 08:45:48 2016
@author: mariarosa
"""

# Import packages
import numpy as np
import csv
from scipy import stats
from sklearn import svm
from sklearn import preprocessing
from sklearn import cross_validation
from sklearn.grid_search import GridSearchCV

# Load entire .csv file but skip header
def load_all_data( data_path):
    f = open(data_path)
    csv_f = csv.reader(f)
    data = []
    for row in csv_f:
        data.append(row)
    data = np.array(data,dtype = float)
    return data

# Load data
data_path = '/Users/maria/Documents/NSPN/rs-fMRI/EVCs_all.csv'
X = load_all_data( data_path)

# Classifier
clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=1./X.shape[1])

# Predictions
preds = []

# LOOCV
print 'LOOCV-UK'
looCV = cross_validation.LeaveOneOut(X.shape[0])
sc = []
for train, test in looCV:
    X_train, X_test = X[train,], X[test,]
    scaler = preprocessing.StandardScaler().fit(X_train)
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)
    clf.fit(X_train)
    preds.append(clf.predict(X_test))
    print clf.predict(X_test)
