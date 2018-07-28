import pandas as pd 
from numpy.linalg.linalg import norm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from os import listdir
from os.path import isfile, join
import sys

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



import matplotlib.pyplot as plt

numslice = int(sys.argv[1])

data_cluster = pd.read_csv('DATA_CLUSTER_WITHOUT_QUERRY.csv')
Data_FOF = pd.read_csv("./FOF_NO_QUERRY_DATA2.csv")

inds = np.loadtxt('index.txt')

#print data_cluster['R'][inds[0,0]],data_cluster['R'][inds[0,1]]

for i in range(12,14):#numslice+1):

	FIG,AX = plt.subplots(1,1,figsize=(8,8))

	R1,R2 = data_cluster['R'][inds[i-1,0]]/1000,data_cluster['R'][inds[i-1,1]]/1000

	id1,id2 = data_cluster['ID'][inds[i-1,0]],data_cluster['ID'][inds[i-1,1]]


	FOF_gal1 = Data_FOF[Data_FOF.grpid ==id1].reset_index()
	FOF_gal2 = Data_FOF[Data_FOF.grpid ==id2].reset_index()
	
	print FOF_gal1.shape,FOF_gal2.shape



	#fig,ax = plt.subplots(1,1,figsize=(8,8))
	fig = plt.figure(figsize=(10,8))
	ax2 = fig.add_subplot(131, projection='3d')
	ax2.set_xlim(0.0,100.0)
	ax2.set_ylim(0.0,100.0)
	ax2.set_zlim(0.0,100.0)

	ax1 = fig.add_subplot(133, projection='3d')
	ax1.set_xlim(0.0,100.0)
	ax1.set_ylim(0.0,100.0)
	ax1.set_zlim(0.0,100.0)

	ax = fig.add_subplot(132, projection='3d')
	ax.set_xlim(0.0,100.0)
	ax.set_ylim(0.0,100.0)
	ax.set_zlim(0.0,100.0)
	



	direc = './SLICED_' + str(i) 
	Params = pd.read_csv(direc+'/Paramters.dat')



	normal = Params['normal']

	R = get_rotation_matrix(normal)
	
	#intersection of the plane
	pplain = (Params['cluster1']+Params['cluster2'])/2.0
	

	Data  = pd.read_csv(direc+'/Final_filament.csv')
	Datap = Transforming(Data,Params)
	P = np.array([Datap['x']-pplain[0],Datap['y']-pplain[1],Datap['z']-pplain[2]])
	Pprime = np.matmul(R,P)

	print "Filament size",Pprime.shape
	
	AX.scatter(Pprime[0,:],Pprime[1,:],s=6,label='Filament',alpha=0.4,color='black')


	Data  = pd.read_csv(direc+'/Final_field.csv')
	Datap = Transforming(Data,Params)
	P = np.array([Datap['x']-pplain[0],Datap['y']-pplain[1],Datap['z']-pplain[2]])
	Pprime = np.matmul(R,P)
	
	AX.scatter(Pprime[0,:],Pprime[1,:],s=6,label='Field',alpha=0.4,color='black')




	
	Data  = pd.read_csv(direc+'/Final_groups.csv')
	Datap = Transforming(Data,Params)
	P = np.array([Datap['x']-pplain[0],Datap['y']-pplain[1],Datap['z']-pplain[2]])
	Pprime = np.matmul(R,P)

	print "Group size",Pprime.shape


	ax.scatter(Datap['x'],Datap['y'],Datap['z'],c='red',s=8)
	ax.scatter(pplain[0],pplain[1],pplain[2],c='black',s=50)
	ax2.scatter(Data['x'],Data['y'],Data['z'],c='blue',s=8)
	ax1.scatter(Pprime[0,:],Pprime[1,:],Pprime[2,:],c='red',s=8)


	#print Pprime[2,:].min(),Pprime[2,:].max()

	AX.scatter(Pprime[0,:],Pprime[1,:],s=6,label='Group',alpha=0.4,color='black')

	C1 = new_point(Params,np.array(Params['cluster1']))
	C2 = new_point(Params,np.array(Params['cluster2']))
	P1 = np.array([C1[0]-pplain[0],C1[1]-pplain[1],C1[2]-pplain[2]])
	P2 = np.array([C2[0]-pplain[0],C2[1]-pplain[1],C2[2]-pplain[2]])
	Cprime1 = np.matmul(R,P1)
	Cprime2 = np.matmul(R,P2)

	

	onlyfiles = [f for f in listdir(direc)]
	

	filament_list = [f for f in onlyfiles if f.startswith('filament_')]

	lengths = np.loadtxt(direc+'/Lengths.dat')
     
	for i,f in enumerate(filament_list):
		dat = np.loadtxt(direc+'/'+f)
		if lengths[i]>1.0:
			#AX.plot(dat[:,0],dat[:,1],color='red')
			AX.plot([],[])

	#AX.scatter(Cprime1[0],Cprime1[1],s=50,color='red')
	#AX.scatter(Cprime2[0],Cprime2[1],s=50,color='red')



	Data_clus = Transforming(FOF_gal1,Params)
	P = np.array([Data_clus['x']-pplain[0],Data_clus['y']-pplain[1],Data_clus['z']-pplain[2]])
	Pprime = np.matmul(R,P)
	#AX.scatter(Pprime[0,:],Pprime[1,:],s=6,alpha=0.4,color='red')

	Data_clus = Transforming(FOF_gal2,Params)
	P = np.array([Data_clus['x']-pplain[0],Data_clus['y']-pplain[1],Data_clus['z']-pplain[2]])
	Pprime = np.matmul(R,P)
	#AX.scatter(Pprime[0,:],Pprime[1,:],s=6,alpha=0.4,color='red')

	AX.set_xlabel('$X^{\prime}$(Mpc)',  fontsize=16)
	AX.set_ylabel('$Y^{\prime}$(Mpc)',  fontsize=16)

	#AX1.legend(loc='upper center',fancybox=True, shadow=True,fontsize=13, ncol=2)
	#AX.legend(fontsize=13)
	AX.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	AX.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
	
	FIG.tight_layout()
	print direc+"/SLICE_colored.png"
	FIG.savefig(direc+"/SLICE__black_colored.png")






