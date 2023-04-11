#!/usr/bin/env python
# coding: utf-8

# # **TPOT AutoML**

# All media | Median molecules log transformed

###############################################################################
## Importing libraries

import os
import time
import numpy as np
from numpy import mean
from numpy import std
from numpy import absolute
import pandas as pd
from joblib import dump, load

from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.metrics import explained_variance_score, max_error, mean_absolute_error, mean_squared_error
from sklearn.metrics import mean_squared_log_error, median_absolute_error, r2_score, accuracy_score
from sklearn.utils import shuffle

from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import make_pipeline, make_union
from sklearn.preprocessing import PolynomialFeatures
from tpot.builtins import StackingEstimator
from tpot.export_utils import set_param_recursive
from sklearn.preprocessing import FunctionTransformer
from copy import copy

from math import sqrt
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt

###############################################################################
## Data loading and transformation

data = pd.read_csv("trainingdata_Yeast_kapp.csv", sep='\t')

col = []
for column in data.columns:
    col.append(column)

target_col = col[4]
features = col[6:len(col)]

X = data[features].values
X[0] = np.log10(abs(X[0]))
X[1] = np.log10(abs(X[1]))

y = data[target_col].values
y = np.log10(y)
y = np.reshape(y, (-1,1))

###############################################################################
## Cross-validation of predictions

model = make_pipeline(
    make_union(
        FunctionTransformer(copy),
        FunctionTransformer(copy)
    ),
    PolynomialFeatures(degree=2, include_bias=False, interaction_only=False),
    StackingEstimator(estimator=KNeighborsRegressor(n_neighbors=49, p=2, weights="distance")),
    ExtraTreesRegressor(bootstrap=True, max_features=0.15000000000000002, min_samples_leaf=1, min_samples_split=9, n_estimators=100)
)

set_param_recursive(model.steps, 'random_state', 137)

X, y = shuffle(X, y)

print("CROSS VALIDATED R2 SCORES: ", cross_val_score(model, X, y.ravel(), scoring = "r2", cv=10, n_jobs = -1, verbose = 1))
