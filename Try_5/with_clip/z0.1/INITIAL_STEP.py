import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import eagleSqlTools as sql
import numpy as np
import pandas as pd
from subprocess import call
from astropy.io import ascii
import os
from numpy.linalg.linalg import norm





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




def RUN_DISPERSE(file,k):
	
	
	call(['/home/ankit/Python_Environments/EAGLE/DisPerSE/bin/delaunay_2D',\
	  file,'-outDir','./SLICED_'+str(k)+'/'])

	
	NDfile = file + '.NDnet'

	call(['/home/ankit/Python_Environments/EAGLE/DisPerSE/bin/mse',\
	  NDfile,'-nsig','3','-upSkl','-forceLoops','-outDir','./SLICED_'+str(k)+'/'])

	SKLfile = NDfile + '_s3.up.NDskl'

	call(['/home/ankit/Python_Environments/EAGLE/DisPerSE/bin/skelconv',\
	  SKLfile,'-breakdown','-assemble','90','-to','vtp','-outDir','./SLICED_'+str(k)+'/'])

	call(['/home/ankit/Python_Environments/EAGLE/DisPerSE/bin/skelconv',\
	  SKLfile,'-to','NDskl_ascii','-outDir','./SLICED_'+str(k)+'/'])


	FILAMENT_FILE = file +'.NDnet_s3.up.NDskl.a.NDskl'

	#np.loadtxt('filament_83.dat')


	with open(FILAMENT_FILE) as infile, open('SLICED_'+str(k)+'/FILAMENTS'+'.dat','w') as outfile:
		copy = False
		for line in infile:
			if line.strip() == "[FILAMENTS]":
				copy = True
			elif line.strip() == "[CRITICAL POINTS DATA]":
				copy = False
			elif copy:
				outfile.write(line)
		    
		#outfile.close()
		#print line

	with open(FILAMENT_FILE) as infile, open('SLICED_'+str(k)+'/Critical_Points' + '.dat', 'w') as outfile:
		copy = False
		for line in infile:
			if line.strip() == "[CRITICAL POINTS]":
				copy = True
			elif line.strip() == "[FILAMENTS]":
				copy = False
			elif copy:
				#print line.split(' ')
				if (not line.startswith(' ')):
					if(len(line)>10):
						outfile.write('%d,%0.6f\n'%(int(line.split(' ')[0]),float(line.split(' ')[1])))




	outfile.close()

	fil_file = open('SLICED_'+str(k)+'/FILAMENTS'+'.dat', 'r')

	lines = fil_file.readlines()

	#print len(lines)

	total_filaments = lines[0]


	l=0
	filRA = []
	filDEC = []
	for j in range(1, len(lines)):
	    #print l
		if lines[j].startswith(' '):
			fp.write(lines[j])

		else :
			l=l+1
			fp = open('SLICED_'+str(k)+'/filament_'+str(l)+'.dat','w')

	fp.close()






def Slice(normal,Data_total,origin_shift,rad,k):


	a,b,c = normal[0],normal[1],normal[2]
	origin_shiftx,origin_shifty,origin_shiftz = origin_shift[0],origin_shift[1],origin_shift[2]


	z_dir=np.array([0.0,0.0,1.0])

	R = get_rotation_matrix(normal)

	trans_matrix = np.zeros((3,Data_total['x'].shape[0]))



	#SHIFT ORIGIN
	trans_matrix[0,:] = Data_total['x']-origin_shiftx
	trans_matrix[1,:] = Data_total['y']-origin_shifty
	trans_matrix[2,:] = Data_total['z']-origin_shiftz



	#print "DISTANCES in P_rot"
	P_rot = np.matmul(R, trans_matrix)


	#print "P_rot shape",P_rot.shape

	#SELECTING THE POINTS
	d = {'x1': P_rot[0, :], 'y1': P_rot[1, :], 'z1': P_rot[2, :],'id':Data_total['id']}
	B = pd.DataFrame(d)

	

	#clipped = np.zeros((B.loc[np.abs(B['z1']) <= rad].shape[1], B.loc[np.abs(B['z1']) <= rad].shape[0]))

	clippedx=[]
	clippedy=[]
	clippedz=[]
	clippedid=[]


	for ii in range(B.shape[0]):
		#print ii
			
		if(np.abs(B['z1'].loc[ii])<=rad):
			clippedx.append(B['x1'].loc[ii])
			clippedy.append(B['y1'].loc[ii])
			clippedz.append(B['z1'].loc[ii])
			clippedid.append(B['id'].loc[ii])


	clipped = np.array([clippedid,clippedx,clippedy,clippedz]).T
        #print clipped
	
	#clipped = clip_data(clipped)

	if (os.path.isdir('./SLICED_'+str(k))):
		print "exits_ SLICED"
	else:
		call(['mkdir','./SLICED_'+str(k)])

	np.savetxt('./SLICED_'+str(k)+'/Clipped.dat',clipped[:,1:3],fmt='%10.4f',header='px py')

        #np.savetxt('./SLICED_'+str(k)+'/Clipped3D.dat',clipped.T)

	Dict = {'id':clipped[:,0],'x':clipped[:,1],'y':clipped[:,2],'z':clipped[:,3]}

	df = pd.DataFrame(Dict)
	df.to_csv('./SLICED_'+str(k)+'/Clipped3D.csv',index=False)

	file = './SLICED_'+str(k)+'/Clipped.dat'
	print "M here"
	RUN_DISPERSE(file,k)



import pandas as pd

Data = ascii.read('DATA_FOR_STRUCT_WITHOUT_QUERRY.csv')
 
data1 = pd.read_csv('CLUSTERS_WITHOUT_QUERRY.csv')

data = np.array([data1['x'].values,data1['y'].values,data1['z'].values]).T


#data = np.loadtxt('Data_Cluster.csv', delimiter=',', skiprows=1)




from itertools import combinations

comb = combinations(range(data.shape[0]), 2)
distlist = []
between1 = []
between2 = []
DMAX = 0.0

for i in comb:
        dist = np.sqrt((data[i[0], 0]-data[i[1], 0])**2 + (data[i[0], 1]-data[i[1], 1])**2 + (data[i[0], 2]-data[i[1], 2])**2 )
                  
       
        if (dist > 0.0):
                distlist.append(dist)
                between1.append(i[0])
                between2.append(i[1])




k=1

idx = np.array([between1,between2]).T

np.savetxt('index.txt',idx,fmt='%d %d')


for i,j in zip(between1,between2):
        print i,j


        cluster1 = np.array([data[i,0],data[i,1],data[i,2]])
        cluster2 = np.array([data[j,0],data[j,1],data[j,2]])


        origin_shift = (cluster1+cluster2)/2.0

        rad = 30.0#10* max(data[i,4],data[j,4])/1000.0

        origin_shiftx, origin_shifty, origin_shiftz = origin_shift[0], origin_shift[1], origin_shift[2]
        #origin_shift = np.array([origin_shiftx,origin_shifty,origin_shiftz])


        cluster_vec = cluster2 - cluster1

        cluster_vec = cluster_vec/norm(cluster_vec)

        norx = np.random.uniform(-1.0,1.0)#np.random.rand()
        nory = np.random.uniform(-1.0,1.0)#np.random.rand()
        norz = -(cluster_vec[0]*norx + cluster_vec[1]*nory)/cluster_vec[2]

        normal = np.array([norx,nory,norz])
        normal = normal/norm(normal)

        
        Parameters = {'normal':normal , 'cluster1':cluster1 , 'cluster2': cluster2}

        D_p = pd.DataFrame(Parameters)


        if os.path.isdir('SLICED_'+str(k)):
        	continue
        else :
        	print "M here"
                call(['mkdir','SLICED_'+str(k)])
                

        D_p.to_csv('SLICED_'+str(k)+'/Paramters.dat',index=False)

        R = get_rotation_matrix(normal)

        #print np.dot(R,cluster2-origin_shift)
       

        #Data_to_slice = Data[:,0:3].T
        Slice(normal,Data,origin_shift,rad,k)

	
        k=k+1

