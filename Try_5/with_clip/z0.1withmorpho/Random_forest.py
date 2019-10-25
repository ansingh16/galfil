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



Total_data = pd.read_csv('/home/ankit/Fortran_code/z0.1/All_Data_Stacked.csv')


Total_data['ENV'][Total_data['d_per'] > 1.0] = 'field'

Total_data['u_minus_r'] = Total_data['u'] - Total_data['r']
Total_data['g_minus_r'] = Total_data['g'] - Total_data['r']

Total_data['sSFR'] = Total_data['SFR']/Total_data['SM']

Total_data['logSM'] = np.log10(Total_data['SM'])

Total_data['State'] = 'Active'

Total_data.loc[Total_data['SFR']*1.0e9/Total_data['SM'] <= 0.01, 'State'] = 'Passive'


features = ['SM','VelDisp','SFR','GM','Tot_Mass',\
			'SM_sh','Metal','Vel','Mass_sh']

Active_gal = Total_data.loc[Total_data['State']=='Active']

X = Active_gal[features].values

y = Active_gal['SFR']



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

'''


# plot feature importance manually
from numpy import loadtxt
from xgboost import XGBClassifier
from matplotlib import pyplot
import numpy as np
from xgboost import plot_importance


def Separate_ENV():

	Total_data = pd.read_csv('/home/ankit/Fortran_code/z0.1/All_Data_Stacked.csv')


	Total_data['ENV'][Total_data['d_per'] > 1.0] = 'field'

	Total_data['u_minus_r'] = Total_data['u'] - Total_data['r']
	Total_data['g_minus_r'] = Total_data['g'] - Total_data['r']

	Total_data['sSFR'] = Total_data['SFR']/Total_data['SM']

	Total_data['logSM'] = np.log10(Total_data['SM'])

	Total_data['State'] = 'Active'

	Total_data.loc[Total_data['SFR']*1.0e9/Total_data['SM'] <= 0.01, 'State'] = 'Passive'


	grouped = Total_data.groupby('ENV')

	for key,data in grouped:
		print key,data.shape
		if key=='filament':
			#data['ENV'][Total_data['d_per'] > 1.0] = 'field'
			Filament = data
		if key=='groups':
			Groups = data
		if key=='field':
			Field = data

	return Groups,Filament,Field



Groups,Filament,Field = Separate_ENV()

print Filament.columns

fig,(ax1) = plt.subplots(1,1,figsize=(6,6))
fig2,(ax2) = plt.subplots(1,1,figsize=(6,6))
fig3,(ax3) = plt.subplots(1,1,figsize=(6,6))


X = Filament[['VelDisp','GM','Tot_Mass','Metal','SM']]
y = Filament['State']
# fit model no training data
model = XGBClassifier()
model.fit(X, y)
# feature importance
print(model.feature_importances_)
# plot
plot_importance(model,ax=ax1)


X = Groups[['VelDisp','GM','Tot_Mass','Metal','SM']]
y = Groups['State']
# fit model no training data
model = XGBClassifier()
model.fit(X, y)
# feature importance
print(model.feature_importances_)
# plot
plot_importance(model,ax=ax2)

X = Field[['VelDisp','GM','Tot_Mass','Metal','SM']]
y = Field['State']
# fit model no training data
model = XGBClassifier()
model.fit(X, y)
# feature importance
print(model.feature_importances_)
# plot
plot_importance(model,ax=ax3)



pyplot.show()
'''
