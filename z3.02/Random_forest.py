import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
from sklearn import metrics
import pandas as pd

from sklearn import preprocessing

from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis as LDA



Total_data = pd.read_csv('All_Data_Stacked.csv',index_col=0)

Total_data['ENV'][Total_data['d_per'] > 1.0] = 'field'

features = ['SM','VelDisp','SFR','GM','Tot_Mass',\
			'SM_sh','Metal','x','y','z','Vel','Mass_sh','u','g','r','i','zmag','Y','J','H','K']



X = Total_data[features].values

y = Total_data['ENV']



from sklearn.cross_validation import cross_val_score, ShuffleSplit
from sklearn.datasets import load_boston
from sklearn.ensemble import RandomForestRegressor
from sklearn import datasets, svm 

clf = svm.SVC(kernel='linear')
clf.fit(X, y)

svm_weights = (clf.coef_ ** 2).sum(axis=0)
svm_weights /= svm_weights.max()

plt.bar(X_indices - .25, svm_weights, width=.2, label='SVM weight',
        color='navy', edgecolor='black')

clf_selected = svm.SVC(kernel='linear')
clf_selected.fit(selector.transform(X), y)

svm_weights_selected = (clf_selected.coef_ ** 2).sum(axis=0)
svm_weights_selected /= svm_weights_selected.max()

plt.bar(X_indices[selector.get_support()] - .05, svm_weights_selected,
        width=.2, label='SVM weights after selection', color='c',
        edgecolor='black')


plt.title("Comparing feature selection")
plt.xlabel('Feature number')
plt.yticks(())
plt.axis('tight')
plt.legend(loc='upper right')
plt.show()