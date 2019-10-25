import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score


import pandas as pd
numslices=55

stack_Data=pd.DataFrame(columns=['id','SM','VelDisp','SFR','GM','Tot_Mass','SM_sh','Metal','x','y','z','Vel',\'Mass_sh','u','g','r','i','zmag','Y','J','H','K','index','ENV','d_long','d_per','d_tot','xslice','yslice'
])

for i in range(1,numslices+1):

    direc = './SLICED_'+ str(i)

    Total_data = pd.read_csv(direc+'/Final_data.csv')

    Total_data['g_minus_r'] = Total_data['g'] - Total_data['r']
    Total_data['u_minus_r'] = Total_data['u'] - Total_data['r']
    Total_data['logSM'] = np.log10(Total_data['SM'])


    stack_Data = stack_Data.append(Total_data)

    

    #Mean_col = Mean_col + Total_data['D'].dropna()


print stack_Data.columns

stack_Data.to_csv('All_Data_Stacked.csv',index=False)
