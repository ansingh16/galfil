import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

from astropy.io import ascii
from astropy.table import Table, join 
import sys

tot = sys.argv[1]

data_all = ascii.read('DATA_FOR_STRUCT_WITHOUT_QUERRY.csv')

for i in range(1,int(tot)+1):

    file1 = 'SLICED_'+str(i)+'/Groups.csv'
    file2 = 'SLICED_'+str(i)+'/Filaments.csv'
    file3 = 'SLICED_'+str(i)+'/Field.csv'

    data_g = ascii.read(file1)
    data_fil = ascii.read(file2)
    data_f = ascii.read(file3)
    
    dat_fil = join(data_fil,data_all,keys='id')
    dat_gr = join(data_g,data_all,keys='id')
    dat_f = join(data_f,data_all,keys='id')

    D2 = dat_fil.to_pandas()
    D3 = dat_f.to_pandas()
    D1 = dat_gr.to_pandas()

    
    
    print D1.shape,D2.shape,D3.shape
 
 
    final_g = 'SLICED_'+str(i)+'/Final_groups.csv'
    final_fil = 'SLICED_'+str(i)+'/Final_filament.csv'
    final_f = 'SLICED_'+str(i)+'/Final_field.csv'



    D1.to_csv(final_g,index=False,float_format='%8.4f')
    D2.to_csv(final_fil,index=False,float_format='%8.4f')
    D3.to_csv(final_f,index=False,float_format='%8.4f')
