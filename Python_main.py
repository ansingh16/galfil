#!/home/ankit/Python_Environments/EAGLE/bin/python
import numpy as np
import Python_part
import sys  
from astropy.table import Table
import pandas as pd
from astropy.io import ascii


tot = int(sys.argv[1])


for i in range(1,tot+1):
  
    print i

    file = './SLICED_'+str(i)+'/Clipped3D.csv'

    data = ascii.read(file)
    
    Python_part.Environment_classifier(data,i)
     
    print "Done Pair: ",i
    '''
    print data_groups.shape
    print data_filament.shape
    print data_field.shape
    


    arr_fil = { 'id':data_filament[:,0].astype(int),\
                'd_per':data_filament[:,3],\
                'd_long':data_filament[:,4],\
                'd_tot':data_filament[:,5]}

    dat = Table(arr_fil)

    data_frame = dat.to_pandas()
    data_frame.set_index('id')



    data_frame.to_csv('./SLICED_'+str(i)+'/Filaments.csv',index=False,float_format='%8.4f')

    arr_Group = {'id':data_groups[:,0].astype(int)}

    dat = Table(arr_Group)

    data_frame = dat.to_pandas()
    data_frame.set_index('id')

    data_frame.to_csv('./SLICED_'+str(i)+'/Groups.csv',index=False)

    arr_field = {'id':data_field[:,0].astype(int)}

    dat = Table(arr_field)

    data_frame = dat.to_pandas()
    data_frame.set_index('id')

    data_frame.to_csv('./SLICED_'+str(i)+'/Field.csv',index=False)

    '''



