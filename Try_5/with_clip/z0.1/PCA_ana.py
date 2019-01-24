from sklearn.preprocessing import StandardScaler

import pandas as pd

Total_data = pd.read_csv('All_Data_Stacked.csv')

features = ['logSM','Metal_1','g_minus_r','u_minus_r','Vel','SFR']

X = Total_data.loc[:, features].values

Y = Total_data.loc[:,['ENV']].values

X = StandardScaler().fit_transform(X)

# Choosing number of Principal Components
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
'''
fig,ax = plt.subplots(1,1,figsize=(6,6)) 
max_pc = 6
pcs = []
totexp_var = []
for i in range(max_pc):
    pca = PCA(n_components=i+1)
    reduced_X = pca.fit_transform(X)
    tot_var = pca.explained_variance_ratio_.sum()
    pcs.append(i+1)
    totexp_var.append(tot_var)
ax.plot(pcs,totexp_var,'r')
ax.plot(pcs,totexp_var,'bs')
ax.set_xlabel('No. of PCs',fontsize = 13)
ax.set_ylabel('Total variance explained',fontsize = 13)
plt.xticks(pcs,fontsize=13)
plt.yticks(fontsize=13)
fig.savefig('PCA_numcomp.png',dpi=600)

'''

from mayavi import mlab

pca = PCA(n_components=3)
reduced_X = pca.fit_transform(X)

print reduced_X.shape[0]
pbar = ProgressBar(widgets=[Percentage(), Bar()],maxval=reduced_X.shape[0]).start()



fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(reduced_X[np.where(Y=='groups'),0],reduced_X[np.where(Y=='groups'),1],c='r',marker='+')

ax.scatter(reduced_X[np.where(Y=='filament'),0],reduced_X[np.where(Y=='filament'),1],c='b',marker='*')

ax.scatter(reduced_X[np.where(Y=='field'),0],reduced_X[np.where(Y=='field'),1],c='green')


ax.set_title('PCA')
plt.show()
