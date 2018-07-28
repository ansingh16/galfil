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

print Total_data.min(),Total_data.max()

# Split the data into 40% test and 60% training
X_train, X_test, y_train, y_test = train_test_split(X, y,test_size=0.30, random_state=123)





std_scale = preprocessing.StandardScaler().fit(X_train)
X_train = std_scale.transform(X_train)
X_test = std_scale.transform(X_test)



lda_clf = LDA()
lda_clf.fit(X_train, y_train)
LDA(n_components=None, priors=None)



pred_train = lda_clf.predict(X_train)

print('Prediction accuracy for the training dataset')
print('{:.2%}'.format(metrics.accuracy_score(y_train, pred_train)))

