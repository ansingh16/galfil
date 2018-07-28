import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np
from os import listdir
from astropy.table import Table, vstack
from astropy.io import ascii

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


def Transforming(Dat,Params):

	
	X1=[]
	Y1=[]
	Z1=[]

	for i in range(Dat.shape[0]):

		r = [Dat[i,0],Dat[i,1],Dat[i,2]]

		rprime = new_point(Params,r)

		X1.append(rprime[0])
		Y1.append(rprime[1])
		Z1.append(rprime[2])


	d = {'x':np.array(X1),'y':np.array(Y1),'z':np.array(Z1)}

	return pd.DataFrame(d)


def Plot_this(DATA,Params,R):
	Dat = np.array([DATA.x,DATA.y,DATA.z]).T
	Datap = Transforming(Dat,Params)
	P = np.array([Datap.x-pplain[0],Datap.y-pplain[1],Datap.z-pplain[2]])
	Pprime = np.matmul(R,P)
	return Pprime


def Mvscol():
			##########################################################M_star vs u-r######################################
			fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(14,6),sharey=True)

			Total_data = pd.read_csv('All_Data_Stacked.csv')
			#Total_data = pd.read_csv('SLICED_13/Final_data.csv')

			Total_data['ENV'][Total_data['ENV']=='field'] = 'tmp1'
			Total_data['ENV'][Total_data['ENV']=='groups'] = 'field'
			Total_data['ENV'][Total_data['ENV']=='tmp1'] = 'groups'

			grouped = Total_data.groupby('ENV')

			for key,data in grouped:
				print key,data.shape

				if key=='filament':
					data['ENV'][Total_data['d_per'] > 1.0] = 'field'
					Filament = data
				if key=='groups':
					Groups = data

				if key=='field':
					Field = data

				
			print Groups.shape,Filament.shape,Field.shape


			cb = ax1.scatter(np.log10(Groups['SM']),Groups['u'] - Groups['r'],c=np.log10(Groups['Metal']/0.012),s=2)
			ax2.scatter(np.log10(Filament['SM']),Filament['u'] - Filament['r'],c=np.log10(Filament['Metal']/0.012),s=2)
			ax3.scatter(np.log10(Field['SM']),Field['u'] - Field['r'],c=np.log10(Field['Metal']/0.012),s=2)
			#cb = ax1.hexbin(np.log10(Groups['SM']),Groups['u'] - Groups['r'],bins='log',gridsize=50,cmap='inferno')
			#ax2.hexbin(np.log10(Groups['SM']),Groups['u'] - Groups['r'],bins='log',gridsize=50,cmap='inferno')
			#ax3.hexbin(np.log10(Groups['SM']),Groups['u'] - Groups['r'],bins='log',gridsize=50,cmap='inferno')



			ax1.set_xticks([8.5,9.5,10.5])
			ax1.set_xlim(8.0,11.5)
			ax1.set_xlabel(r'$\mathrm{log M_{*}}$',fontsize=16)
			ax1.set_ylabel(r'$\mathrm{u - r}$',fontsize=16)
			ax1.set_title('Groups',fontsize=16)
			ax1.minorticks_on()
			ax1.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
			ax1.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
			ax2.set_xticks([8.5,9.5,10.5])
			ax2.set_xlim(8.0,11.5)
			ax2.set_xlabel(r'$\mathrm{log M_{*}}$',fontsize=16)
			ax2.set_title('Filament',fontsize=16)
			ax2.minorticks_on()
			ax2.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
			ax2.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
			ax3.set_xticks([8.5,9.5,10.5])
			ax3.set_xlim(8.0,11.5)
			ax3.set_xlabel(r'$\mathrm{log M_{*}}$',fontsize=16)
			ax3.set_title('Field',fontsize=16)
			ax3.minorticks_on()
			ax3.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
			ax3.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
			ax = [ax1,ax2,ax3]
			fig.subplots_adjust(wspace=0.0, hspace=0.0)
			C = fig.colorbar(cb, ax=ax)
			C.set_label(r'$Z/Z_{\odot}$',fontsize=14)
			C.ax.tick_params(labelsize=14) 

			#counts,xbins,ybins=np.histogram2d(np.log10(Groups['SM']),Groups['u'] - Groups['r'],bins=50)
			# make the contour plot
			'''
			mylevels=[100.0, 1000.0]

			ax3.contour(counts.transpose(),mylevels,extent=[xbins.min(),xbins.max(),
			    ybins.min(),ybins.max()],linewidths=2,colors='black',
			    linestyles='solid')
			ax2.contour(counts.transpose(),mylevels,extent=[xbins.min(),xbins.max(),
			    ybins.min(),ybins.max()],linewidths=2,colors='black',
			    linestyles='solid')
			ax1.contour(counts.transpose(),mylevels,extent=[xbins.min(),xbins.max(),
			    ybins.min(),ybins.max()],linewidths=2,colors='black',
			    linestyles='solid')
			'''
			fig2,ax2 = plt.subplots()

			ax2.scatter(Groups['xslice'],Groups['yslice'],s=2,color='red')
			ax2.scatter(Filament['xslice'],Filament['yslice'],s=2,color='blue')
			ax2.scatter(Field['xslice'],Field['yslice'],s=2,color='green')

			fig.savefig('slice13_Stellar_mass_vs_uminusr.png',dpi=600,bbox_inches='tight')

			plt.show()

def historgram():
			######################################################HISTOGRAM PLOT#############################################


			Total_data = pd.read_csv('All_Data_Stacked.csv',index_col=0)
			Total_data['ENV'][Total_data['d_per'] > 1.0] = 'field'
			fig2,(ax2) = plt.subplots(1,1,figsize=(6,6))
			Groups = Total_data[Total_data['ENV']=='groups']
			Filament = Total_data[Total_data['ENV']=='filament']
			Field = Total_data[Total_data['ENV']=='field']



			y1,binEdges1=np.histogram(Groups['u'] - Groups['r'],bins='auto',density=True)
			y2,binEdges2=np.histogram(Filament['u'] - Filament['r'],bins='auto',density=True)
			y3,binEdges3=np.histogram(Field['u'] - Field['r'],bins='auto',density=True)

			bincenters1 = 0.5*(binEdges1[1:]+binEdges1[:-1])
			bincenters2 = 0.5*(binEdges2[1:]+binEdges2[:-1])
			bincenters3 = 0.5*(binEdges3[1:]+binEdges3[:-1])

			ax2.plot(binEdges1[:-1], y1,'-o', color='m', linewidth=2, ms=8, label='Groups')
			ax2.plot(binEdges2[:-1], y2,'-o', color='peru', linewidth=2, ms=8, label='Filament')
			ax2.plot(binEdges3[:-1], y3,'-o', color='r', linewidth=2, ms=8   ,label='Field')  

			ax2.set_xlabel(r'$\mathrm{u - r}$')
			ax2.set_ylabel('Probability density')

			fig2.legend()
			plt.show()

def jointplot():

			####################################################JOINT PLOT#####################################################



			fig = plt.figure(figsize=(8, 8))
			grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
			main_ax = fig.add_subplot(grid[:-1, 1:])
			y_hist = fig.add_subplot(grid[:-1, 0], xticklabels=[], sharey=main_ax)
			x_hist = fig.add_subplot(grid[-1, 1:], yticklabels=[], sharex=main_ax)
			main_ax.plot(np.log10(Filament['SM']),Filament['u'] - Filament['r'],'ok', markersize=3, alpha=0.2)
			x_hist.hist(np.log10(Filament['SM']), 40, histtype='stepfilled',orientation='vertical', color='gray')
			x_hist.invert_yaxis()
			y_hist.hist(Filament['u'] - Filament['r'], 40, histtype='stepfilled',orientation='horizontal', color='gray')
			y_hist.invert_xaxis()
			plt.show()



def Voilin_plot():
			##############################################VIOLIN PLOT#############################################################


			import seaborn as sns
			sns.set(style='white',font_scale=1.5)
			fig,ax = plt.subplots(1,1,figsize=(8,8))
			Total_data = pd.read_csv('All_Data_Stacked.csv')
			#Total_data = pd.read_csv('SLICED_13/Final_data.csv')

			Total_data['ENV'][Total_data['ENV']=='field'] = 'tmp1'
			Total_data['ENV'][Total_data['ENV']=='groups'] = 'field'
			Total_data['ENV'][Total_data['ENV']=='tmp1'] = 'groups'

			Total_data['u_minus_r'] = Total_data['u'] - Total_data['r']

			ax = sns.violinplot(x="ENV", y="u_minus_r",
			                     data=Total_data, palette="muted", split=True)

			ax.set(ylabel=r'$\mathrm{u - r}$',xlabel='Environment')
			fig.savefig('violinplot.png',dpi=600)
			plt.show()


def Colored_slice():

			################################################################################################################



			FIG1,AX1 = plt.subplots(1,1,figsize=(8,8))
			global direc
			direc = './SLICED_' + str(13) 
			'''
			Params = pd.read_csv(direc+'/Paramters.dat')
			normal = Params['normal']
			global pplain 
			pplain = (Params['cluster1']+Params['cluster2'])/2.0

			R = get_rotation_matrix(normal)
				

			Total_data = pd.read_csv(direc+'/Final_data.csv')
			#Total_data['ENV'][Total_data['d_per'] > 1.0] = 'field'
			for ind,row in Total_data.iterrows():

				if(row['d_per']>1.0):
					Total_data.loc[ind, 'ENV'] = 'field' 

			'''
	
			dat = pd.read_csv(direc+'/All_env.csv')

			
			
			grouped = dat.groupby('ENV')

			for key,data in grouped:
				print key,data.shape
				
				if key=='filament':
					data['ENV'][data['d_per'] > 1.0] = 'field'
					Filament = data
				if key=='group':
					Groups = data
				if key=='field':
					Field = data
			
			AX1.scatter(Field['xslice'],Field['yslice'],s=2,color='green',label='field')
			AX1.scatter(Filament['xslice'],Filament['yslice'],s=2,color='blue',label='filament')
			AX1.scatter(Groups['xslice'],Groups['yslice'],s=2,color='red',label='group')

			'''

			for ind,row in data.iterrows():

				if(row['d_per']>1.0):
					data.loc[ind, 'ENV'] = 'field' 
			
			
			
			onlyfiles = [f for f in listdir(direc)]
				

			filament_list = [f for f in onlyfiles if f.startswith('filament_')]

			lengths = np.loadtxt(direc+'/Lengths.dat')
			 
			for q,f in enumerate(filament_list):
					dat = np.loadtxt(direc+'/'+f)
					if lengths[q]>=1.0:
						if dat.shape[0]>3:
							AX1.plot(dat[:,0],dat[:,1],color='red')
							AX1.plot([],[])

			'''
			
			
			AX1.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
			AX1.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)


			AX1.set_xlabel(r'$\mathrm{X^{\prime}}$',fontsize=16)
			AX1.set_ylabel(r'$\mathrm{Y^{\prime}}$',fontsize=16)
			AX1.legend(fontsize=13)
			FIG1.savefig(direc+'/Colored_without_filament.png',dpi=600)
			

			'''
			onlyfiles = [f for f in listdir(direc)]
				

			filament_list = [f for f in onlyfiles if f.startswith('filament_')]

			lengths = np.loadtxt(direc+'/Lengths.dat')
			     
			for q,f in enumerate(filament_list):
					dat = np.loadtxt(direc+'/'+f)
					if lengths[q]>1.0:
						AX1.plot(dat[:,0],dat[:,1],color='red')
						AX1.plot([],[])

						



			FIG1.legend(fontsize=13)

			FIG1.tight_layout()

			plt.show()
		'''


def All_data():
	D_final = Table()

	for i in range(1,56):

		direc = 'SLICED_'+str(i)

		Di = ascii.read(direc+'/Final_data.csv',names=['id','SM','VelDisp','SFR','GM',\
			'Tot_Mass','SM_sh','Metal','x','y','z','Vel','Mass_sh','u','g','r','i',\
			'zmag','Y','J','H','K','index','ENV','d_long','d_per','d_tot','xslice','yslice'])
		
		D_final = vstack([D_final, Di])

		print len(D_final['d_per'])


	#ascii.write('All_Data_Stacked.csv',D_final)
	Write_DF = D_final.to_pandas()

	Write_DF.to_csv('All_Data_Stacked.csv',index=False)


def dvscol():

	fig,ax = plt.subplots(1,1,figsize=(8,8))

	Total_data = pd.read_csv('All_Data_Stacked.csv')

	Total_data = Total_data[Total_data['ENV']=='filament']

	Total_data.sort_values('d_per',inplace=True)

	Total_data = Total_data.reset_index(drop=True)

	Total_data['g_minus_r'] = Total_data['g'] - Total_data['r']

	feature=['d_per','g_minus_r']

	stack_Data = Total_data[feature]


	stack_Data = stack_Data.dropna()
	stack_Data.sort_values('d_per',inplace=True)

	stack_Data.to_csv('Stacked_Data.csv',index=False)

	stack_Data['D'] =  stack_Data.rolling(6000).median()['d_per']
	stack_Data['col'] =  stack_Data.rolling(6000).median()['g_minus_r']


	stack_Data.plot(ax = ax,x='D',y='col',legend=False,linewidth=1.5,color='black')

	ax.set_xlabel(r'$\mathrm{d_{per}}$(Mpc)',fontsize=16)
	ax.set_ylabel(r'$\mathrm{g - r}$',fontsize=16)
	ax.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	ax.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)

	#ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

	fig.savefig('d_per_vs_gminusr.png',dpi=600)
	plt.show()


def main():

	print "Please enter an option:\n"
	print "\n 0: Get stacked file\n"
	print "\n 1: M_star vs u-r colored by metallicity\n"
	print "\n 2: histogram of u-r for different env\n"
	print "\n 3: jointplot for one env \n"
	print "\n 4. Voilin_plot for u-r different Environment\n"
	print "\n 5. Colored_slice for one slice\n"
	print "\n 6. d vs color plot\n"

	opt=input('option?')

	if opt==0:
		All_data()
	if opt==1:
		Mvscol()
	if opt==2:
		histogram()
	if opt==3:
		jointplot()
	if opt==4:
		Voilin_plot()
	if opt==5:
		Colored_slice()
	if opt==6:
		dvscol()


if __name__=='__main__':

	main()