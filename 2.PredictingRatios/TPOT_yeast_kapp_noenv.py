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
from sklearn.metrics import explained_variance_score, max_error, mean_absolute_error, mean_squared_error
from sklearn.metrics import mean_squared_log_error, median_absolute_error, r2_score, accuracy_score

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

data = pd.read_csv("trainingdata_Yeast_kapp_noenv.csv", sep='\t')

col = []
for column in data.columns:
    col.append(column)

target_col = col[4]
features = col[5:len(col)]

X = data[features].values
#X[0] = np.log10(abs(X[0]))
#X[1] = np.log10(abs(X[1]))

y = data[target_col].values
y = np.log10(y)
y = np.reshape(y, (-1,1))

###############################################################################
## Model configuration

model = make_pipeline(
    make_union(
        FunctionTransformer(copy),
        FunctionTransformer(copy)
    ),
    PolynomialFeatures(degree=2, include_bias=False, interaction_only=False),
    StackingEstimator(estimator=KNeighborsRegressor(n_neighbors=49, p=2, weights="distance")),
    ExtraTreesRegressor(bootstrap=True, max_features=0.15000000000000002, min_samples_leaf=1, min_samples_split=9, n_estimators=100)
)
# Fix random state for all the steps in exported pipeline
set_param_recursive(model.steps, 'random_state', 137)

###############################################################################
## Fit model

X_train, X_valid, y_train, y_valid = train_test_split(X, y, test_size=0.2, random_state=137)

model.fit(X_train, y_train.ravel())

###############################################################################
## Model evaluation

y_rescaled = y_valid

predict_valid = model.predict(X_valid)
predict_valid = np.reshape(predict_valid, (-1,1))

#baseline_preds = y_rescaled[:,target_col.index("ratio")]
#baseline_errors = abs(baseline_preds - y_rescaled)
#errors = abs(predict_valid - y_rescaled)
#mape = 100 * (errors / y_rescaled)
#accuracy = 100 - np.mean(mape)


print('\n')
print("MODEL VALIDATION METRICS", "\n")

#print("Average baseline error: ", round(np.mean(baseline_errors),2))
#print("Mean absolute error: ", round(np.mean(errors),2))
#print("Accuracy: ", round(accuracy, 2), "%", "\n")

print("Explained variance regression score: ", explained_variance_score(y_rescaled, predict_valid))
print("R2 score: ", r2_score(y_rescaled, predict_valid))
print("Adjusted R2 score: ", 1 - (1-r2_score(y_rescaled, predict_valid))*(len(y)-1)/(len(y)-X.shape[1]-1), "\n")

print("Maximum residual error: ", max_error(y_rescaled, predict_valid))
print("Median absolute error: ", median_absolute_error(y_rescaled, predict_valid))
print("Mean absolute error: ", mean_absolute_error(y_rescaled, predict_valid))
print("Mean squared error: ", mean_squared_error(y_rescaled, predict_valid))
print("Root mean squared error:", sqrt(mean_squared_error(y_rescaled, predict_valid)))
#print("Mean squared logarithmic error: ", mean_squared_log_error(y_rescaled, predict_valid))

###############################################################################
## Correlation between experimental data and predicted values

pearson = stats.pearsonr(y_rescaled.ravel(), predict_valid.ravel())
spearman = stats.spearmanr(y_rescaled.ravel(), predict_valid.ravel())

print('\n')
print('CORRELATION BETWEEN EXPERIMENTAL DATA AND PREDICTED VALUES', '\n')
print('Pearson\'s r:', pearson[0], 'p-value:', pearson[1])
print('Spearman\'s r:', spearman[0], 'p-value:', spearman[1], '\n')


plot_data = pd.DataFrame()
plot_data['Known abundance'] = y_rescaled.ravel()
plot_data['Predicted abundance'] = predict_valid.ravel()

sns.regplot(x='Known abundance', y='Predicted abundance', data=plot_data)

plt.savefig('TPOT_Yeast_kapp_noenv_valid_correlation.jpg', format='jpg', dpi=120)


residuals = y_rescaled - predict_valid

plot_data = pd.DataFrame()
plot_data['Residuals'] = residuals.ravel()
plot_data['Predicted abundance'] = predict_valid.ravel()

sns.jointplot(y='Residuals', x='Predicted abundance', data=plot_data, kind="resid")

plt.savefig('TPOT_Yeast_kapp_noenv_valid_residuals.jpg', format='jpg', dpi=120)

###############################################################################
## Predicted values

#predict_valid = np.expm1(predict_valid)
#y_rescaled = np.expm1(y_rescaled)

predict_valid = 10**predict_valid
y_rescaled = 10**y_rescaled

fmt = '%-8s%-20s%s'

print('VALIDATION DATASET PROTEIN ABUNDANCE PREDICTIONS', '\n')

print(fmt % ('', 'Eval data', 'Prediction'))
for i, (eval_row, pred_row) in enumerate(zip(y_rescaled, predict_valid)):
    print(fmt % (i, eval_row, pred_row))


###############################################################################
## Model saving
filename = 'TPOT_Yeast_kapp_noenv_exported_pipeline.sav'
dump(model, filename)

