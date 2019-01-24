import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import numpy as np 
import sys
from astropy.io import ascii
from astropy.table import Table, join
import pandas as pd
import matplotlib.pyplot as plt
from os import listdir


numslice = int(sys.argv[1])

def get_rotation_matrix(i_v, unit=None):
	
	# From http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q38
	if unit is None:
		unit = [0.0, 0.0, 1.0]
	# Normalize vector length
	i_v /= np.linalg.norm(i_v)
	# Get axis
	uvw = np.cross(i_v, unit)

	
	# compute trig values - no need to go through arccos and back
	rcos = np.dot(i_v, unit)
	rsin = np.linalg.norm(uvw)
    #normalize and unpack axis
	if not np.isclose(rsin, 0):
		uvw /= rsin
	
	
	u,v,w = uvw

	# Compute rotation matrix - re-expressed to show structure
	return (
        rcos * np.eye(3) +
        rsin * np.array([
            [ 0, -w,  v],
            [ w,  0, -u],
            [-v,  u,  0]
        ]) +
        (1.0 - rcos) * uvw[:,None] * uvw[None,:]
	)


def new_point(Params,r):

	normal = Params['normal']
	
	pplain = (Params['cluster1']+Params['cluster2'])/2.0

	v = r - pplain

	dist = np.dot(v,normal)


	rprime = r - dist*normal #r - np.dot((r - Params['cluster1']),normal)*normal

	#rprime = r_0 + t_1*e_1 + t_2*e_2 + s*normal

	return rprime


def Transforming(Data,Params):

	
	X1=[]
	Y1=[]
	Z1=[]

	for i in range(Data.shape[0]):

		r = [Data['x'][i],Data['y'][i],Data['z'][i]]

		rprime = new_point(Params,r)

		X1.append(rprime[0])
		Y1.append(rprime[1])
		Z1.append(rprime[2])


	d = {'x':np.array(X1),'y':np.array(Y1),'z':np.array(Z1)}

	return pd.DataFrame(d)



def Plot_this(DATA,Params,R):
	Datap = Transforming(DATA,Params)
	P = np.array([Datap['x']-pplain[0],Datap['y']-pplain[1],Datap['z']-pplain[2]])
	Pprime = np.matmul(R,P)
	return Pprime


Available_data = ascii.read('DATA_FOR_ANALYSIS_WITHOUT_QUERRY.csv',names=['id','SM','VelDisp','SFR',\
	'GM','Tot_Mass','SM_sh','Metal','SF_Metal','NSF_Metal',\
	'SF_O','SF_H','NSF_O','NSF_H','Stars_O','Stars_H','x','y','z','SubGrpNum','Vel','Mass_sh','grpid','u','g','r','i','zmag','Y','J','H','K','M_dot'])
print Available_data.keys()



for i in range(1,numslice+1):

	FIG1,AX1 = plt.subplots(1,1,figsize=(8,8))
	FIG2,AX2 = plt.subplots(1,1,figsize=(8,8))



	file_g = './SLICED_' + str(i) + '/Final_groups.csv'
	Data_g = pd.read_csv(file_g)
	Data_g['ENV'] = 'groups'

	file_fil = './SLICED_' + str(i) + '/Final_filament.csv'
	Data_fil = pd.read_csv(file_fil)
	Data_fil['ENV'] = 'filament'

	file_f = './SLICED_' + str(i) + '/Final_field.csv'
	Data_f = pd.read_csv(file_f)
	Data_f['ENV'] = 'field'

	f = [Data_g,Data_f,Data_fil]
	
	Data_all_env = pd.concat(f).reset_index(drop=True)



	Data_table = Table.from_pandas(Data_all_env)

	All = Available_data.to_pandas()
	All = All.round({'x': 4, 'y': 4,'z': 4,'Metal':4})


	Available_data = Table.from_pandas(All)

	#print Available_data.keys()

	T = join(Available_data,Data_table,keys=['id','x','y','z','SM','Metal'])

	DF = T.to_pandas()

	DF.to_csv('./SLICED_'+str(i)+'/Final_data.csv',index=False)

	print "Total galaxies:",DF.shape

	'''

	#PLOT THE GALAXIES FOR WHICH THE DATA IS AVAILABLE

	direc = './SLICED_' + str(i) 
	Params = pd.read_csv(direc+'/Paramters.dat')
	normal = Params['normal']
	pplain = (Params['cluster1']+Params['cluster2'])/2.0

	R = get_rotation_matrix(normal)
	
	Pprime= Plot_this(Data_all_env,Params,R)
	#AX1.scatter(Pprime[0,:],Pprime[1,:],s=4,alpha=0.4,color='black')
	AX2.scatter(Pprime[0,:],Pprime[1,:],s=4,alpha=0.4,color='black')


	G= DF[DF['ENV']=='groups'].reset_index()
	Fil= DF[DF['ENV']=='filament'].reset_index()
	F= DF[DF['ENV']=='field'].reset_index()



	Pprime = Plot_this(Fil,Params,R)

	
	print "Filament galaxies:",Pprime.shape

	AX1.scatter(Pprime[0,:],Pprime[1,:],s=6,alpha=0.4,color='blue',label='Filament Galaxies')
	Pprime = Plot_this(F,Params,R)

	print "Field galaxies:",Pprime.shape


	AX1.scatter(Pprime[0,:],Pprime[1,:],s=6,alpha=0.4,color='green',label='Field Galaxies')
	Pprime = Plot_this(G,Params,R)
	AX1.scatter(Pprime[0,:],Pprime[1,:],s=6,alpha=0.4,color='red',label='Group')

	print "Group galaxies:",Pprime.shape

	
	onlyfiles = [f for f in listdir(direc)]
	

	filament_list = [f for f in onlyfiles if f.startswith('filament_')]

	lengths = np.loadtxt(direc+'/Lengths.dat')
     
	for q,f in enumerate(filament_list):
		dat = np.loadtxt(direc+'/'+f)
		if lengths[q]>1.0:
			AX1.plot(dat[:,0],dat[:,1],color='red')
			AX1.plot([],[])

			AX2.plot(dat[:,0],dat[:,1],color='red')
			AX2.plot([],[])
	


	#Data_test_groups = pd.read_csv('./GROUPS_VERIFY_WITHOUT_QUERRY.csv')
	Data_test_groups = pd.read_csv('./GROUPS_STRUCT_WITHOUT_QUERRY.csv')
	Data_test_groups = Data_test_groups.round({'x': 4, 'y': 4,'z': 4,'Metal':4})
	Data_slice = pd.read_csv('./SLICED_'+str(i)+'/Clipped3D.csv')	
	Data_slice['id'] = Data_slice['id'].astype(int)
	Dat_G = join(Table.from_pandas(DF),Table.from_pandas(Data_test_groups),keys=['id','x','y','z'])
	Group = Dat_G.to_pandas()

	Pprime = Plot_this(Group,Params,R)
	#AX2.scatter(Pprime[0,:],Pprime[1,:],s=6,alpha=1.0,color='red',label='Groups')






	#PLOT GROUPS
	#AX1.set_title('Galaxies for which data is available')
	AX1.set_xlabel('$X^{\prime}$(Mpc)',  fontsize=16)
	AX1.set_ylabel('$Y(Mpc)^{\prime}$',  fontsize=16)

	#AX1.legend(loc='upper center',fancybox=True, shadow=True,fontsize=13, ncol=2)
	AX1.legend(fontsize=13)
	AX1.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	AX1.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
	#FIG1.tight_layout()
	FIG1.savefig(direc+'/SLICED_'+str(i)+'_Magnitude_available.png',dpi=600)

	#AX2.set_title('Galaxies classified as groups by FOF')
	AX2.set_xlabel('$X^{\prime}$(Mpc)',  fontsize=16)
	AX2.set_ylabel('$Y^{\prime}$(Mpc)',  fontsize=16)
	AX2.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	AX2.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)

	AX2.legend(fontsize=13)
	FIG1.tight_layout()

	FIG2.savefig(direc+'/SLICED_'+str(i)+'_FOF_Groups.png',dpi=600)
	'''



