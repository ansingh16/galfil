import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Data_FOF = pd.read_csv("./FOF_NO_QUERRY_DATA2.csv")

Data_FOF['logSM'] = np.log10(Data_FOF['astar'])

#ax = plt.subplots(3,4,figsize=(14,20))
hist = Data_FOF.hist(column='logSM',by='grpid',normed=True,figsize=(10,10),sharey=True,sharex=True)


for i in range(len(hist)):
    axes = hist[i]

    for ax in axes:
        ax.set_xlabel(r'$log M_{*} M_{\odot}$')
        ax.set_ylabel('Freq')
        

plt.subplots_adjust(wspace=0)
plt.savefig('Cluster_histogram_amass.png')
plt.show()
