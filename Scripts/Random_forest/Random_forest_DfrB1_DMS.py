#!/usr/bin/env python
# coding: utf-8

# # Random forest DMS
# 
# This script runs the random forest model on the data from the differences in fitness effects: deltaS_weak (S_weak - S_opt)


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge
from sklearn.tree import DecisionTreeRegressor
from IPython.display import HTML

from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

from sklearn.preprocessing import PolynomialFeatures
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from statistics import mean

from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import mean_squared_error
from sklearn.metrics import accuracy_score

from rfpimp import permutation_importances
from sklearn.base import clone 

from sklearn import tree
import graphviz
from sklearn.tree import _tree
import sklearn
sklearn.__version__


dataset = pd.read_csv('../../Data/Complete_datasets/dataset_diffNorm_ara0.2_ara0.01_index_differences.txt', sep='\t')
dataset

# Remove stop codons and rows that are not present in the PDB structure
dataset = dataset[(dataset['Residue'] != '*') & 
                  (pd.notna(dataset['Mean_ddG_stab_HET'])) &
                  (pd.notna(dataset['diffNormScore']))
                 ]

df = dataset.drop(['Entropy','Position','WT_Residue', 'Residue', 'Arabinose',
                   'Secondary_structure','Solvent_accessibility', 'mean_sel_coeff', 
                   'mean_sel_coeff_2', 'Arabinose_2'],axis=1)

df


X = df.drop(['diffNormScore'],axis=1)
Y = df['diffNormScore']

# Normalize all the features to the range of -1 , 1
X = X / X.max()


## Train a random forest model
X_train,X_test,y_train,y_test =train_test_split(X,Y,test_size=0.2, random_state=100)

model_rf = RandomForestRegressor(n_estimators=500, oob_score=True, random_state=100)
model_rf.fit(X_train, y_train) 
pred_train_rf= model_rf.predict(X_train)

print('Mean squared error (train):', np.sqrt(mean_squared_error(y_train,pred_train_rf)))
print('R2 score (train):', r2_score(y_train, pred_train_rf))

pred_test_rf = model_rf.predict(X_test)

print('Mean squared error (test):', np.sqrt(mean_squared_error(y_test,pred_test_rf)))
print('R2 score (test):', r2_score(y_test, pred_test_rf))

## Train the random forest again but adding a random variable
np.random.seed(100)
X['random_var'] = np.random.normal(loc = 0, scale = 1, size = X.shape[0])

X_train,X_test,y_train,y_test =train_test_split(X,Y,test_size=0.2, random_state = 100)

model_rf = RandomForestRegressor(n_estimators=500, oob_score=True, random_state=100)
model_rf.fit(X_train, y_train) 
pred_train_rf= model_rf.predict(X_train)

print('Mean squared error (train):', np.sqrt(mean_squared_error(y_train,pred_train_rf)))
print('R2 score (train):', r2_score(y_train, pred_train_rf))

pred_test_rf = model_rf.predict(X_test)

print('Mean squared error (test):', np.sqrt(mean_squared_error(y_test,pred_test_rf)))
print('R2 score (test):', r2_score(y_test, pred_test_rf))


# Use cross-validation on this preliminary model
cross_val_n = 5
print('Five-fold cross validation of the random forest model:')
cross_validations = cross_val_score(estimator = model_rf, X = X_train, y = y_train, cv=cross_val_n , scoring = r2)
print(cross_validations)
print(np.mean(cross_validations), np.std(cross_validations) / np.sqrt(cross_val_n))
print('------')

cross_validations = cross_val_score(estimator = model_rf, X = X_train, y = y_train, cv=5, scoring = r2)
print(cross_validations)

## Check cross-validation accuracy
cross_validations_pred = cross_val_predict(estimator = model_rf, X = X_train, y = y_train, cv=5)

get_ipython().run_line_magic('matplotlib', 'inline')
## Scatterplot of the random forest predictions 
fig = plt.figure()

ax = fig.add_axes([0,0,1,1])
ax.scatter(cross_validations_pred, y_train)
ax.set_xlabel('Predicted deltaS (random forest, cross-validations)', fontsize = 20)
ax.set_ylabel('Observed deltaS', fontsize = 20)
cross_validations_pred[1, ]
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(14)

plt.show()


# Figures for the accuracy of the predictions and selecting the best variables
get_ipython().run_line_magic('matplotlib', 'inline')
## Scatterplot of the random forest predictions 
fig = plt.figure()

ax = fig.add_axes([0,0,1,1])
ax.scatter(pred_test_rf, y_test)
ax.set_xlabel('Predicted fitness effects (random forest)', fontsize = 20)
ax.set_ylabel('Observed fitness effects', fontsize = 20)

for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(14)

plt.show()

## Save the corresponding files for predictions on validation set and test set
## Results with all variables
df_pred_test = pd.DataFrame(list(zip(y_test,  pred_test_rf)), columns = ['test_data', 'pred_data'])
df_pred_test.to_csv('../../Data/Random_forest_results/diffNorm_ara0.2_ara0.01/pred_rf_allVariables.txt', sep = '\t')

## Results of predictions in the cross-validation
df_pred_crossval = pd.DataFrame(list(zip(y_train, cross_validations_pred)), columns = ['test_data', 'pred_data'])
df_pred_crossval.to_csv('../../Data/Random_forest_results/diffNorm_ara0.2_ara0.01/crossval_rf_allVariables.txt', sep = '\t')


# ## Feature selection
# Define a function to use permutation to estimate relative importances
def r2(rf, X_train, y_train):
    return r2_score(y_train, rf.predict(X_train))

# Use permutation to estimate relative importances
perm_imp_rfpimp = permutation_importances(model_rf, X_train, y_train, r2)

get_ipython().run_line_magic('matplotlib', 'inline')
fig = plt.figure()

ax = fig.add_axes([0,0,1,3])
ax.barh(list(perm_imp_rfpimp.index), perm_imp_rfpimp['Importance'])
ax.set_xlabel('Relative importance', fontsize = 20)
ax.set_ylabel('Feature', fontsize = 20)

for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(14)

plt.show()

best_features = perm_imp_rfpimp[perm_imp_rfpimp['Importance'] >= perm_imp_rfpimp['Importance']['random_var']]
best_features

new_X = X[list(best_features.index)]

# Train a new random forest with the selected variables and the random variable
X_train,X_test,y_train,y_test =train_test_split(new_X,Y,test_size=0.2, random_state = 100)

model_rf = RandomForestRegressor(n_estimators=500, oob_score=True, random_state=100)
model_rf.fit(X_train, y_train) 
pred_train_rf= model_rf.predict(X_train)

print('Mean squared error (train):', np.sqrt(mean_squared_error(y_train,pred_train_rf)))
print('R2 score (train):', r2_score(y_train, pred_train_rf))

pred_test_rf = model_rf.predict(X_test)

print('Mean squared error (test):', np.sqrt(mean_squared_error(y_test,pred_test_rf)))
print('R2 score (test):', r2_score(y_test, pred_test_rf))

get_ipython().run_line_magic('matplotlib', 'inline')
## Scatterplot of the random forest predictions 
fig = plt.figure()

ax = fig.add_axes([0,0,1,1])
ax.scatter(pred_test_rf, y_test)
ax.set_xlabel('Predicted fitness effects (random forest)', fontsize = 20)
ax.set_ylabel('Observed fitness effects', fontsize = 20)

for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(14)

plt.show()

# Cross-validation
cross_val_n = 5

print('Five-fold cross validation of the random forest model:')
cross_validations = cross_val_score(estimator = model_rf, X = X_train, y = y_train, cv=cross_val_n, scoring = r2)
print(cross_validations)
print(np.mean(cross_validations), np.std(cross_validations) / np.sqrt(cross_val_n))
print('------')

cross_validations_pred = cross_val_predict(estimator = model_rf, X = X_train, y = y_train, cv=cross_val_n)

get_ipython().run_line_magic('matplotlib', 'inline')
## Scatterplot of the random forest predictions 
fig = plt.figure()

ax = fig.add_axes([0,0,1,1])
ax.scatter(cross_validations_pred, y_train)
ax.set_xlabel('Predicted deltaS (random forest, cross-validations)', fontsize = 20)
ax.set_ylabel('Observed deltaS', fontsize = 20)
cross_validations_pred[1, ]
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(14)

plt.show()


# Since this is a simpler model, we can test relative importance by leaving one
# variable out and retraining
# Function from https://explained.ai/rf-importance/index.html#intro
def dropcol_importances(rf, X_train, y_train):
    rf_ = clone(rf)
    rf_.random_state = 100
    rf_.fit(X_train, y_train)
    baseline = rf_.oob_score_
    imp = []
    for col in X_train.columns:
        X = X_train.drop(col, axis=1)
        rf_ = clone(rf)
        rf_.random_state = 100
        rf_.fit(X, y_train)
        o = rf_.oob_score_
        imp.append(baseline - o)
    imp = np.array(imp)
    I = pd.DataFrame(
            data={'Feature':X_train.columns,
                  'Importance':imp})
    I = I.set_index('Feature')
    I = I.sort_values('Importance', ascending=True)
    return I

importances_drop_col = dropcol_importances(model_rf, X_train, y_train)

get_ipython().run_line_magic('matplotlib', 'inline')
fig = plt.figure()

ax = fig.add_axes([0,0,1,1])
ax.barh(list(importances_drop_col.index), importances_drop_col['Importance'])
ax.set_xlabel('Relative importance', fontsize = 16)
ax.set_ylabel('Feature', fontsize = 16)

for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(12)

plt.show()

### Save the tables
## Relative importances permutation (all)
perm_imp_rfpimp.to_csv('../../Data/Random_forest_results/model_diffFit_permImportances_allVariables.txt', sep = '\t')

## Predictions for test set(best variables)
df_pred_test_best = pd.DataFrame(list(zip(y_test,  pred_test_rf)), columns = ['test_data', 'pred_data'])
df_pred_test_best.to_csv('../../Data/Random_forest_results/pred_rf_bestVariables.txt', sep = '\t')

## Predictions for cross-validation (best variables)
df_pred_crossval_best = pd.DataFrame(list(zip(y_train, cross_validations_pred)), columns = ['test_data', 'pred_data'])
df_pred_crossval_best.to_csv('../../Data/Random_forest_results/crossval_rf_bestVariables.txt', sep = '\t')

## Relative importances drop column (best variables)
importances_drop_col.to_csv('../../Data/Random_forest_results/model_diffFit_dropCol_bestVariables.txt', sep = '\t')

