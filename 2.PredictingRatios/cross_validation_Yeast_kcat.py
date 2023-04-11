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

from sklearn.linear_model import RidgeCV, SGDRegressor
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline, make_union
from sklearn.preprocessing import Normalizer
from tpot.builtins import StackingEstimator
from xgboost import XGBRegressor
from tpot.export_utils import set_param_recursive

from math import sqrt
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt

###############################################################################
## Data loading and transformation

data = pd.read_csv("trainingdata_Yeast_kcat.csv", sep='\t')

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
    StackingEstimator(estimator=RidgeCV()),
    Normalizer(norm="l1"),
    StackingEstimator(estimator=SGDRegressor(alpha=0.001, eta0=1.0, fit_intercept=True, l1_ratio=0.0, learning_rate="invscaling", loss="huber", penalty="elasticnet", power_t=100.0)),
    XGBRegressor(learning_rate=0.1, max_depth=7, min_child_weight=15, n_estimators=100, n_jobs=1, objective="reg:squarederror", subsample=0.8, verbosity=0)
)

set_param_recursive(model.steps, 'random_state', 137)

X, y = shuffle(X, y)

print("CROSS VALIDATED R2 SCORES: ", cross_val_score(model, X, y.ravel(), scoring = "r2", cv=10, n_jobs = -1, verbose = 1))
