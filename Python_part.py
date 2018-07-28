import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.constants import c
from astropy.cosmology import WMAP9 as cosmo
from astropy.coordinates import search_around_3d
from numpy.linalg import norm
from sklearn.cluster import DBSCAN
from sklearn import metrics
import pandas as pd
from scipy.spatial import distance
#from sympy import Point, Line
import sympy
from shapely.geometry import LineString,Point
from shapely.ops import nearest_points
from shapely.geometry import MultiPoint
from astropy.io import ascii
from astropy.table import Table, join
import collections
import matplotlib.pyplot as plt
from os import listdir

radius = 3.0
limit_log = 3.0

def closest_node(node, nodes):
	
	
	closest_index = distance.cdist([node], nodes).argmin()
	return nodes[closest_index]



def clip_data(data):
	data = data[data[:, 3] > min(data[:, 3] + 5.0)]
	data = data[data[:, 3] < max(data[:, 3] - 5.0)]

	data = data[data[:, 4] > min(data[:, 4] + 5.0)]
	data = data[data[:, 4] < max(data[:, 4] - 5.0)]
	return data


def clustering( r, sam):

	dbsc = DBSCAN(eps = r, min_samples = sam).fit(file_RADEC_array)
	labels = dbsc.labels_
	core_samples = np.zeros_like(labels, dtype = bool)
	core_samples[dbsc.core_sample_indices_] = True
	labels = dbsc.labels_
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

	return core_samples,labels,n_clusters_


def Environment_classifier(Data, K):
	
	
	
	direc = 'SLICED_'+str(K)
	fig,AX1 = plt.subplots(1,1,figsize=(8,8))

	points = np.loadtxt(direc+'/Clipped3D.csv',skiprows=1,delimiter=',')

	all_data = points[:,1:]

	radius=3.0
	#data = np.random.rand(100,3)
	#filament_data = np.random.rand(100,3)

	onlyfiles = [f for f in listdir(direc)]
	filament_list = [f for f in onlyfiles if f.startswith('filament_')]

	file = direc+'/Clipped3D.csv'

	Data = ascii.read(file)
	    
	Groups_data = pd.read_csv('./GROUPS_STRUCT_WITHOUT_QUERRY.csv')

	Tot = Data['id'].tolist()
	Tot = [int(i) for i in Tot]


	model_groups = collections.Counter(Groups_data.grpid)
	new_df = pd.concat( [Groups_data.query('grpid==%d'%key) for key,val in model_groups.items() if val > 4 ] ) 


	G = new_df['id'].tolist()

	labels=[]
	Glabel=[]

	for q in range(len(Data['x'])):
	    if(Tot[q] not in G):
	            labels.append(q)
	    else:
	            Glabel.append(q)

	print len(Glabel)

	data = all_data[labels,:]	#Data1[labels,:]
	IDnogrp = points[labels,0]
	IDgrp = points[Glabel,0]



	lengths = np.loadtxt(direc+'/Lengths.dat')
	Sky_all_points = SkyCoord(data[:, 0],data[:, 1],0.0,unit='Mpc', representation='cartesian')
	all_filament_gal=[]
	all_filament_gal_dist=[]
	all_filament_gal_dist_fil=[]
	all_filament_gal_dist_tot_fil=[]


	for q,f in enumerate(filament_list):
	    dat = np.loadtxt(direc+'/'+f)
	    filament_data=np.zeros((dat.shape[0],3))
	    filament_data[:,0] = dat[:,0]
	    filament_data[:,1] = dat[:,1]
	    
	    if lengths[q]>=1.0:
	        if dat.shape[0]>3:
	            Sky_all_points = SkyCoord(data[:, 0],data[:, 1],0.0,unit='Mpc', representation='cartesian')

	            Filament_all_points = SkyCoord(filament_data[:, 0],filament_data[:, 1],0.0,unit='Mpc',representation='cartesian')
	            fil_value, sky_values, sep, dist = search_around_3d(Filament_all_points, Sky_all_points,distlimit=radius * u.Mpc)

	            dat1 = np.array([fil_value, sky_values, sep, dist]).T
	            dat1 = dat1[dat1[:, 3].argsort()]
	            data_fil_ind = dat1[:, 0].astype(int)
	            data_gal_ind = dat1[:, 1].astype(int)
	            data_dist_gal = dat1[:, 3]
	            nearby_gal, indexes = np.unique(data_gal_ind, return_index=True)
	            near_fil_point = data_fil_ind[indexes]
	            dist_fil = data_dist_gal[indexes]
	            
	            d_per = np.zeros((len(nearby_gal), 1))
	            d_clus = np.zeros((len(nearby_gal), 1))
	            d_totlen = np.zeros((len(nearby_gal), 1))

	            x_fil = filament_data[:, 0]
	            y_fil = filament_data[:, 1]

	            list1 = []

	            for m in range(filament_data.shape[0]):
	                list1.append(Point(x_fil[m], y_fil[m]))

	            for j in range(len(nearby_gal)):
	                p = Point(data[nearby_gal[j], 0], data[nearby_gal[j], 1])

	                p1 = (data[nearby_gal[j], 0], data[nearby_gal[j], 1])

	                D = radius
	                
	                L = LineString(list1)

	                NP = L.interpolate(L.project(p))

	                route = MultiPoint(list1)
	                
	                NV = nearest_points(route, p)[0]

	                indp = list1.index((NV))

	                dist1=  np.sqrt((NV.x - NP.x)**2+(NV.y - NP.y)**2)
	                fildist = dist1
	                if(indp<(len(list1)-indp)):
	                        for m in range(0, indp):
	                            x_along_fil1 = x_fil[m]
	                            y_along_fil1 = y_fil[m]

	                            x_along_fil2 = x_fil[m + 1]
	                            y_along_fil2 = y_fil[m + 1]

	                            dist1 = np.sqrt((x_along_fil1 - x_along_fil2)**2+(x_along_fil1 - x_along_fil2)**2)

	                            fildist += dist1

	                elif(indp > (len(list1)-indp)):
	                        for m in range(indp, len(list1) - 1):
	                            x_along_fil1 = x_fil[m]
	                            y_along_fil1 = y_fil[m]

	                            x_along_fil2 = x_fil[m + 1]
	                            y_along_fil2 = y_fil[m + 1]
	                            dist1 = np.sqrt((x_along_fil1 - x_along_fil2)**2+(x_along_fil1 - x_along_fil2)**2)

	                            fildist += dist1

	                elif (indp == len(list1) / 2):

	                        for m in range(indp, len(list1) - 1):
	                            x_along_fil1 = x_fil[m]
	                            y_along_fil1 = y_fil[m]

	                            x_along_fil2 = x_fil[m + 1]
	                            y_along_fil2 = y_fil[m + 1]

	                            dist1 = np.sqrt((x_along_fil1 - x_along_fil2)**2+(x_along_fil1 - x_along_fil2)**2)

	                            fildist += dist1

	                
	                
	                x_gal = float(sympy.N(p.x))
	                y_gal = float(sympy.N(p.y))

	                dist = np.sqrt((x_gal - NP.x)**2+(y_gal - NP.y)**2)
	                d_per[j, 0] = dist
	                
	                d_clus[j, 0] = fildist

	                d_totlen[j, 0] = lengths[q - 1]

	                
	                
	                if (nearby_gal[j] in all_filament_gal):
	                        ind = all_filament_gal.index(nearby_gal[j])
	                        if (all_filament_gal_dist[ind] < dist):
	                            continue
	                        else:
	                            all_filament_gal_dist[ind] = dist
	                            all_filament_gal_dist_fil[ind] = fildist
	                            all_filament_gal_dist_tot_fil[ind] = lengths[q - 1]
	                            all_filament_gal[ind] = nearby_gal[j]
	                else:

	                        all_filament_gal.append(nearby_gal[j])
	                        all_filament_gal_dist.append(dist)
	                        all_filament_gal_dist_fil.append(fildist)
	                        all_filament_gal_dist_tot_fil.append(lengths[q - 1])

	            
	            AX1.plot(dat[:,0],dat[:,1],color='red')
	            AX1.plot([],[])

	data_field = np.array([data[i,:] for i in range(data.shape[0]) if i not in all_filament_gal])
	id_f = np.array([IDnogrp[i] for i in range(data.shape[0]) if i not in all_filament_gal])
	id_fil = IDnogrp[all_filament_gal]

	id_f = [int(w) for w in id_f.tolist()]
	id_fil = [int(w) for w in id_fil.tolist()]
	IDgrp = [int(w) for w in IDgrp.tolist()]

	

	data_groups = all_data[Glabel,:]


	dist_array = np.array([all_filament_gal_dist,all_filament_gal_dist_fil,all_filament_gal_dist_tot_fil]).T
	
	data_filament = data[all_filament_gal,:]

	data_filament = np.append(data[all_filament_gal,:],dist_array,axis=1)

	dictf = {'id':id_f,'xslice':data_field[:,0],'yslice':data_field[:,1]}
	dictG = {'id':IDgrp,'xslice':data_groups[:,0],'yslice':data_groups[:,1]}
	dictfil = {'id':id_fil,'xslice':data_filament[:,0],'yslice':data_filament[:,1],\
	            'd_per':np.array(all_filament_gal_dist),'d_long':np.array(all_filament_gal_dist_fil),\
	            'd_tot':np.array(all_filament_gal_dist_tot_fil)}



	DG = pd.DataFrame(dictG)
	Df = pd.DataFrame(dictf)
	Dfil = pd.DataFrame(dictfil)
	
	Dfil.set_index('id')
	Df.set_index('id')
	DG.set_index('id')

	Dfil.to_csv(direc+'/Filaments.csv',index=False)
	DG.to_csv(direc+'/Groups.csv',index=False)
	Df.to_csv(direc+'/Field.csv',index=False)
	


	Dfil['ENV'] = 'filament'
	Df['ENV'] = 'field'
	DG['ENV'] = 'group'
	D_total = pd.concat([DG,Df,Dfil])
	D_total.to_csv(direc+'/All_env.csv',index=False)

	

	print DG.shape,Dfil.shape,Df.shape
	AX1.scatter(data[all_filament_gal,0],data[all_filament_gal,1],color='blue',s=2)
	AX1.scatter(all_data[Glabel,0],all_data[Glabel,1],s=2,color='red')
	AX1.scatter(Df['xslice'],Df['yslice'],s=2,color='green')
	AX1.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	AX1.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
	AX1.set_xlabel(r'$\mathrm{X^{\prime}}$',fontsize=16)
	AX1.set_ylabel(r'$\mathrm{Y^{\prime}}$',fontsize=16)

	fig.savefig(direc+'/Colored.png',dpi=600)
	







#plt.show()




