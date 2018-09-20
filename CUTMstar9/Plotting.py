import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
from scipy import stats

import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np
from os import listdir
from astropy.table import Table, vstack
from astropy.io import ascii
import ConfigParser
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata

#plt.style.use('dark_background')

Path = 'Figures/'

Colors = ['#ED5752','#4D648D','#A1BE95']

#Colors = ['white','white','white']

Config = ConfigParser.ConfigParser()

Config.read("./Config_plotting.ini")


def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1


def binning(col, cut_points, labels=None):
		  #Define min and max values:
		  minval = col.min()
		  maxval = col.max()
		  #minval = 8.0#col.min()
		  #maxval = 11.0#col.max()


		  #create list by adding min and max to cut_points
		  break_points = [minval] + cut_points + [maxval]

		  #if no labels provided, use default labels 0 ... (n-1)
		  if not labels:
		    labels = range(len(cut_points)+1)

		  #Binning using cut function of pandas
		  colBin = pd.cut(col,bins=break_points,labels=labels,include_lowest=True)

		  return colBin

		

def grid(x, y, z, resX=5000, resY=5000):
    "Convert 3 column data to matplotlib grid"
    xi = linspace(min(x), max(x), resX)
    yi = linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi,interp='linear')
    X, Y = meshgrid(xi, yi)
    return X, Y, Z

def Separate_ENV():

	Total_data = pd.read_csv('All_Data_Stacked.csv')
	
	
	Total_data['ENV'][Total_data['d_per'] > 1.0] = 'field'

	Total_data['u_minus_r'] = Total_data['u'] - Total_data['r']
	Total_data['g_minus_r'] = Total_data['g'] - Total_data['r']

	Total_data['sSFR'] = Total_data['SFR']/Total_data['SM']

	Total_data['logSM'] = np.log10(Total_data['SM'])

	Total_data['State'] = 'Active'

	Total_data.loc[Total_data['SFR']*1.0e9/Total_data['SM'] <= 0.01, 'State'] = 'Passive'


	grouped = Total_data.groupby('ENV')

	for key,data in grouped:
		print key,data.shape
		if key=='filament':
			#data['ENV'][Total_data['d_per'] > 1.0] = 'field'
			Filament = data
		if key=='groups':
			Groups = data
		if key=='field':
			Field = data

	return Groups,Filament,Field


def Preprocess(Groups1,Filament1,Field1,xprop1,yprop1,remove0x,remove0y,logx,logy):

	if remove0x=='y':
				
				Groups1 = Groups1.loc[Groups1[xprop1] != 0.0 ]
				Filament1 = Filament1.loc[Filament1[xprop1] != 0.0 ]
				Field1 = Field1.loc[Field1[xprop1] != 0.0]
			
	if remove0y=='y':
				print "M HERE"
				Groups1 = Groups1.loc[Groups1[yprop1] != 0.0 ]
				Filament1 = Filament1.loc[Filament1[yprop1] != 0.0 ]
				Field1 = Field1.loc[Field1[yprop1] != 0.0]
			
	print "This is minimum",Groups1[yprop1].min()
	
	if logx=='y':
				Groups1['log'+xprop1] = np.log10(Groups1[xprop1])
				Filament1['log'+xprop1] = np.log10(Filament1[xprop1])
				Field1['log'+xprop1] = np.log10(Field1[xprop1])


	if logy=='y':
				Groups1['log'+yprop1] = np.log10(Groups1[yprop1])
				Filament1['log'+yprop1] = np.log10(Filament1[yprop1])
				Field1['log'+yprop1] = np.log10(Field1[yprop1])

	
	
	return Groups1,Filament1,Field1


def _plot_settings_panel(AX,title=' ',xlabel=' ',ylabel=' '):

	AX.set_xlabel(xlabel,fontsize=16)
	AX.set_ylabel(ylabel,fontsize=16)
	AX.set_title(title,fontsize=16)
	AX.minorticks_on()
	AX.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	AX.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
	AX.label_outer()

def Mvscol():
			##########################################################M_star vs u-r######################################
			fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(14,6),sharey=True)

			Groups,Filament,Field = Separate_ENV()

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

			fig.savefig(Path+'Stellar_mass_vs_uminusr.png',dpi=600,bbox_inches='tight')

			plt.show()

def histogram(prop,logit,dens):
			######################################################HISTOGRAM PLOT#############################################

			fig2,ax2 = plt.subplots(1,1,figsize=(8,8))

			Groups,Filament,Field = Separate_ENV()
			
			if logit=='y':
				Groups = Groups[Groups[prop]!=0.0]
				Filament = Filament[Filament[prop]!=0.0]
				Field = Field[Field[prop]!=0.0]

				Groups['log'+prop] = np.log10(Groups[prop])
				Filament['log'+prop] = np.log10(Filament[prop])
				Field['log'+prop] = np.log10(Field[prop])
				prop = 'log'+prop

			if dens=='y':
				dens1 = True
			else:
				dens1 = False


			y1,binEdges1=np.histogram(Groups[prop],bins='auto',density=dens1)
			y2,binEdges2=np.histogram(Filament[prop],bins='auto',density=dens1)
			y3,binEdges3=np.histogram(Field[prop],bins='auto',density=dens1)

			bincenters1 = 0.5*(binEdges1[1:]+binEdges1[:-1])
			bincenters2 = 0.5*(binEdges2[1:]+binEdges2[:-1])
			bincenters3 = 0.5*(binEdges3[1:]+binEdges3[:-1])

			ax2.plot(binEdges1[:-1], y1,'-', color='m', linewidth=2, ms=8, label='Groups')
			ax2.plot(binEdges2[:-1], y2,'--', color='peru', linewidth=2, ms=8, label='Filament')
			ax2.plot(binEdges3[:-1], y3,'-.', color='r', linewidth=2, ms=8   ,label='Field')  

			ax2.set_xlabel(prop,fontsize=14)
			ax2.set_ylabel('Probability density',fontsize=14)
			ax2.tick_params(axis='both', which='minor', length=3, width=2, labelsize=12)
			ax2.tick_params(axis='both', which='major', length=5, width=2, labelsize=12)
			

			ax2.legend(fontsize=14)
			plt.show()

def jointplot(props,xlabel=None,ylabel=None,remove0x=None,remove0y=None,logx=None,logy=None):

			xprop,yprop = props


			####################################################JOINT PLOT#####################################################

			Colors = ['#3F6C45','#CB6318','#29a8ab']

			Groups_tmp,Filament_tmp,Field_tmp = Separate_ENV()

			Groups,Filament,Field = Preprocess(Groups_tmp,Filament_tmp,Field_tmp,\
				xprop,yprop,remove0x,remove0y,logx,logy)
			
			Groups = Groups.sample(1000)
			Filament = Filament.sample(1000)
			Field = Field.sample(1000)


			if logx=='y':
				xprop = 'log'+xprop
			if logy=='y':
				yprop = 'log'+yprop

			fig = plt.figure(figsize=(8, 8))
			grid = plt.GridSpec(4, 4, hspace=0.0, wspace=0.0)
			main_ax = fig.add_subplot(grid[:-1, 1:])
			y_hist = fig.add_subplot(grid[:-1, 0], sharey=main_ax)
			x_hist = fig.add_subplot(grid[-1, 1:], sharex=main_ax)
			
			
			
			main_ax.plot(Groups[xprop],Groups[yprop],'o',color=Colors[2], markersize=3, alpha=1.0,label='Group')
			main_ax.plot(Filament[xprop],Filament[yprop],'o',color=Colors[1], markersize=3, alpha=1.0,label='Filament')
			main_ax.plot(Field[xprop],Field[yprop],'o',color=Colors[0], markersize=3, alpha=1.0,label='Field')
			
			main_ax.get_xaxis().set_visible(False)
			main_ax.get_yaxis().set_visible(False)

			
			y_hist.set_ylabel(r'$\mathrm{'+ylabel+'}$',fontsize=14)

			x_hist.set_xlabel(r'$\mathrm{'+xlabel+'}$',fontsize=14)
			

			y_hist.tick_params(axis='both', which='minor', length=3, width=2, labelsize=12)
			y_hist.tick_params(axis='both', which='major', length=5, width=2, labelsize=12)
			
			x_hist.tick_params(axis='both', which='major', length=5, width=2, labelsize=12)
			x_hist.tick_params(axis='both', which='minor', length=3, width=2, labelsize=12)
			
			#y_hist.set_ylim(0.01,)
			#x_hist.set_ylim(0.01,)


			x_hist.hist([Field[xprop],Filament[xprop],Groups[xprop]],\
				color=Colors,histtype='bar',bins=10,orientation='vertical',align='mid',normed=True)
			

			x_hist.invert_yaxis()
			y_hist.hist([Field[yprop],Filament[yprop],Groups[yprop]],\
				color=Colors,histtype='bar',bins=10,orientation='horizontal',align='mid',normed=True)
			y_hist.invert_xaxis()


			y_hist.xaxis.set_ticks_position('top')
			'''
			y_hist.set_yticks(pos_binsy)
			x_hist.set_xticks(pos_binsx)

			'''
			leg = main_ax.legend(fontsize=14)

			for lh in leg.legendHandles: 
				lh._legmarker.set_alpha(1.0)
				lh._legmarker.set_markersize(5)


			fig.savefig(Path+xprop+'vs'+yprop+'.png',dpi=600)
			
	
			plt.show()
			
			
def Voilin_plot():
			##############################################VIOLIN PLOT#############################################################
			Total_data = pd.read_csv('All_Data_Stacked.csv')
			Total_data['u_minus_r'] = Total_data['u'] - Total_data['r']


			import seaborn as sns
			sns.set(font_scale=1.5)
			sns.set_style('white')

			fig,ax = plt.subplots(1,1,figsize=(8,8))

			Groups,Filament,Field = Separate_ENV()
			ax = sns.violinplot(x="ENV", y="u_minus_r",
			        data=Total_data, palette="muted", split=True)

			ax.set(ylabel=r'$\mathrm{u - r}$',xlabel='Environment')
			fig.savefig(Path+'violinplot.png',dpi=600)
			plt.show()


def Colored_slice(sliceno):

			################################################################################################################



			FIG1,AX1 = plt.subplots(1,1,figsize=(8,8))
			FIG2,AX2 = plt.subplots(1,1,figsize=(8,8))
			FIG3,AX3 = plt.subplots(1,1,figsize=(8,8))


			direc = './SLICED_' + str(sliceno) 

			'''
		
			global direc
			
			Params = pd.read_csv(direc+'/Paramters.dat')
			normal = Params['normal']
			global pplain 
			pplain = (Params['cluster1']+Params['cluster2'])/2.0

			R = get_rotation_matrix(normal)
			'''
				

			Total_data = pd.read_csv(direc+'/Final_data.csv')
			Total_data['ENV'][Total_data['d_per'] > 1.0] = 'field'



			grouped = Total_data.groupby('ENV')

			for key,data in grouped:
				print key,data.shape
				if key=='filament':
					#data['ENV'][Total_data['d_per'] > 1.0] = 'field'
					Filament = data
				if key=='groups':
					Groups = data
				if key=='field':
					Field = data


	

			#Groups,Filament,Field = Preprocess(Groups,Filament,Field,xprop,yprop,'y','y',logx,logy)


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
			
			'''

			AX1.scatter(Field['xslice'],Field['yslice'],marker='s',s=10,facecolor='none',edgecolor='k',label='field')
			AX1.scatter(Filament['xslice'],Filament['yslice'],marker='^',s=10,color='r',label='filament')
			AX1.scatter(Groups['xslice'],Groups['yslice'],s=6,facecolor='black',edgecolor='k',label='group',alpha=0.3)

			
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
			
			AX2.scatter(Field['xslice'],Field['yslice'],s=4,color='grey',label='field',alpha=0.1)
			AX2.scatter(Filament['xslice'],Filament['yslice'],s=4,color='grey',label='filament',alpha=0.1)
			AX2.scatter(Groups['xslice'],Groups['yslice'],s=4,color='grey',label='group',alpha=0.1)

			AX3.scatter(Field['xslice'],Field['yslice'],s=4,color='k',label='field',alpha=0.5)
			AX3.scatter(Filament['xslice'],Filament['yslice'],s=4,color='k',label='filament',alpha=0.5)
			AX3.scatter(Groups['xslice'],Groups['yslice'],s=4,color='k',label='group',alpha=0.5)


			'''
			for ind,row in data.iterrows():

				if(row['d_per']>1.0):
					data.loc[ind, 'ENV'] = 'field' 
			
			
			'''
			
			onlyfiles = [f for f in listdir(direc)]
				

			filament_list = [f for f in onlyfiles if f.startswith('filament_')]

			lengths = np.loadtxt(direc+'/Lengths.dat')
			 
			for q,f in enumerate(filament_list):
					dat = np.loadtxt(direc+'/'+f)
					if lengths[q]>=1.0:
						if dat.shape[0]>3:
							AX2.plot(dat[:,0],dat[:,1],color='red')
							AX2.plot([],[])
			

			
			AX1.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
			AX1.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
			AX1.set_xlabel(r'$\mathrm{X^{\prime}}$',fontsize=16)
			AX1.set_ylabel(r'$\mathrm{Y^{\prime}}$',fontsize=16)
			

			AX2.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
			AX2.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
			AX2.set_xlabel(r'$\mathrm{X^{\prime}}$',fontsize=16)
			AX2.set_ylabel(r'$\mathrm{Y^{\prime}}$',fontsize=16)

			AX3.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
			AX3.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
			AX3.set_xlabel(r'$\mathrm{X^{\prime}}$',fontsize=16)
			AX3.set_ylabel(r'$\mathrm{Y^{\prime}}$',fontsize=16)
			
			

			FIG1.legend(fontsize=14,markerscale=2)

		
			FIG1.savefig(Path+'/Color_without_filament.png',dpi=600)
			FIG2.savefig(Path+'/No_color_with_filament.png',dpi=600)
			FIG3.savefig(Path+'/No_color_without_filament.png',dpi=600)

			plt.show()
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
			'Tot_Mass','SM_sh','Metal','SF_Metal','NSF_Metal','SF_O','SF_H','x','y','z','SubGrpNum','Vel','Mass_sh','grpid','u','g','r','i',\
			'zmag','Y','J','H','K','index','ENV','d_long','d_per','d_tot','xslice','yslice'])
		
		D_final = vstack([D_final, Di])

		print len(D_final['d_per'])


	#ascii.write('All_Data_Stacked.csv',D_final)
	Write_DF = D_final.to_pandas()

	Write_DF.to_csv('All_Data_Stacked.csv',index=False)


def dvscol():

	fig,ax = plt.subplots(1,1,figsize=(6,6))

	Total_data = pd.read_csv('All_Data_Stacked.csv')

	#Total_data.drop_duplicates(subset=['id','ENV'],inplace=True)


	Total_data = Total_data[Total_data['ENV']=='filament']

	Total_data.sort_values('d_per',inplace=True)

	Total_data = Total_data.reset_index(drop=True)

	Total_data['g_minus_r'] = Total_data['g'] - Total_data['r']

	print Total_data['d_per'].max()


	feature=['d_per','g_minus_r']

	stack_Data = Total_data[feature]


	stack_Data = stack_Data.dropna()
	#stack_Data.sort_values('d_per',inplace=True)

	#stack_Data.to_csv('Stacked_Data.csv',index=False)

	stack_Data['D'] =  stack_Data.rolling(150).median()['d_per']
	stack_Data['col'] =  stack_Data.rolling(150).median()['g_minus_r']


	stack_Data.plot(ax = ax,x='D',y='col',legend=False,linewidth=1.5,color='k')

	ax.set_xlabel(r'$\mathrm{d_{per}}$(Mpc)',fontsize=16)
	ax.set_ylabel(r'$\mathrm{g - r}$',fontsize=16)
	ax.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	ax.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)

	#ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
	fig.tight_layout()
	fig.savefig('d_per_vs_gminusr.png',dpi=600)
	plt.show()

def Scatter(props,colprop=None,xlabel=None,ylabel=None,logx=None,logy=None,logcol=None):

	fig,ax = plt.subplots(1,1,figsize=(8,8))

	xprop,yprop = props

	Total_data = pd.read_csv('All_Data_Stacked.csv')
	
	Total_data['g_minus_r'] = Total_data['g'] - Total_data['r']
	Total_data['u_minus_r'] = Total_data['u'] - Total_data['r']

	
	Total_data['ENV'][Total_data['d_per'] > 1.0] = 'field'


	if logx=='y':
				Total_data = Total_data.loc[Total_data[xprop] != 0.0 ]
				Total_data['log'+xprop] = np.log10(Total_data[xprop])
				xprop = 'log'+xprop


	if logy=='y':
				Total_data = Total_data.loc[Total_data[yprop] != 0.0 ]
			
				Total_data['log'+yprop] = np.log10(Total_data[yprop])
				yprop = 'log'+yprop

	if colprop is not None:

				if logcol=='y':

					Total_data = Total_data.loc[Total_data[xprop] != 0.0 ]

					Total_data['log'+colprop] = np.log10(Total_data[colprop])
					colprop = 'log'+colprop

				ax.scatter(Total_data[xprop],Total_data[yprop],c=Total_data[colprop],s=2,alpha=0.2)
	else:

				ax.scatter(Total_data[xprop],Total_data[yprop],color='black',s=2,alpha=0.2)


	ax.set_xlabel(r'$\mathrm{'+xlabel+'}$',fontsize=14)
	ax.set_ylabel(r'$\mathrm{'+ylabel+'}$',fontsize=14)


	ax.tick_params(axis='both', which='minor', length=3, width=2, labelsize=12)
	ax.tick_params(axis='both', which='major', length=5, width=2, labelsize=12)
			
	ax.tick_params(axis='both', which='major', length=5, width=2, labelsize=12)
	ax.tick_params(axis='both', which='minor', length=3, width=2, labelsize=12)
			
	plt.show()	


def Contour(props,xlabel=None,ylabel=None,logx=None,logy=None):
	
	
	#import matplotlib.patches as mpatches
	from matplotlib.lines import Line2D

	import seaborn as sns; sns.set(color_codes=True)
	sns.set_context(context='paper')
	sns.set_style("white")

	#plt.style.use("dark_background")

	#sns.set_style("white")

	new_plt = sns.light_palette("green",as_cmap=True)
	plt_fil = sns.dark_palette("purple", as_cmap=True)


	xprop,yprop = props

	Groups,Filament,Field = Separate_ENV()

	Groups,Filament,Field = Preprocess(Groups,Filament,Field,xprop,yprop,'y','y',logx,logy)


	g = sns.JointGrid(x='log'+xprop, y='log'+yprop, data=Groups,height=6)
	
	sns.kdeplot(
	    Filament['log'+xprop], Filament['log'+yprop], cmap="Greys",
	    shade=False, shade_lowest=False, ax=g.ax_joint,
	    linestyles='--',label='Field',n_levels=5
	)
	
	sns.kdeplot(
	    Groups['log'+xprop], Groups['log'+yprop], cmap="Greys",
	    shade=False, shade_lowest=False, ax=g.ax_joint,
	    linestyles='-',label='Group',n_levels=5
	)
	

	sns.distplot(
	    Groups['log'+xprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls ='-'), ax=g.ax_marg_x,axlabel=False
	)
	sns.distplot(
	    Filament['log'+xprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls='--'), ax=g.ax_marg_x,axlabel=False
	)
	sns.distplot(
	    Groups['log'+yprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls ='-'), ax=g.ax_marg_y, vertical=True,axlabel=False
	)
	sns.distplot(
	    Filament['log'+yprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls='--'), ax=g.ax_marg_y, vertical=True,axlabel=False
	)
	
	#legend_properties = {'size':12,'color':['b','r']}
	#legendMain=g.ax_joint.legend(prop=legend_properties,loc='upper right')

	#patch1 = mpatches.Patch(color='r', label='Group')
	#patch2 = mpatches.Patch(color='b', label='Filament')

	#all_handles = (patch1, patch2)

	#leg = g.ax_joint.legend(handles=all_handles)
	#g.ax_joint.add_artist(leg)

	g.ax_joint.set_xlabel(r"$\mathrm{"+xlabel+"}$",fontsize=12)
	g.ax_joint.set_ylabel(r"$\mathrm{"+ylabel+"}$",fontsize=12)

	g.ax_joint.minorticks_on()
	g.ax_joint.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	g.ax_joint.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
	
	

	######################################################################

	f = sns.JointGrid(x='log'+xprop, y='log'+yprop, data=Filament,height=6)
	
	sns.kdeplot(
	    Filament['log'+xprop], Filament['log'+yprop], cmap="Greys",
	    shade=False, shade_lowest=False, ax=f.ax_joint,
	    linestyles='--',label='Filament',n_levels=5
	)
	sns.kdeplot(
	    Field['log'+xprop], Field['log'+yprop], cmap="Greys",
	    shade=False, shade_lowest=False, ax=f.ax_joint,
	    linestyles='-',label='Field',n_levels=5
	)
	sns.distplot(
	    Filament['log'+xprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls ='--'), ax=f.ax_marg_x,axlabel=False
	)
	sns.distplot(
	    Field['log'+xprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls='-'), ax=f.ax_marg_x,axlabel=False
	)
	sns.distplot(
	    Filament['log'+yprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls ='--'), ax=f.ax_marg_y, vertical=True,axlabel=False
	)
	sns.distplot(
	    Field['log'+yprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls='-'), ax=f.ax_marg_y, vertical=True,axlabel=False
	)
	
	#legend_properties = {'size':12,'color':['b','r']}
	#legendMain=f.ax_joint.legend(prop=legend_properties,loc='upper right')

	#patch2 = mpatches.Patch(color='b', label='Filament')
	#patch1 = mpatches.Patch(color='r', label='Field')

	#all_handles = (patch1, patch2)

	#leg = f.ax_joint.legend(handles=all_handles)
	#f.ax_joint.add_artist(leg)

	f.ax_joint.set_xlabel(r"$\mathrm{"+xlabel+"}$",fontsize=12)
	f.ax_joint.set_ylabel(r"$\mathrm{"+ylabel+"}$",fontsize=12)
	
	f.ax_joint.minorticks_on()
	f.ax_joint.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	f.ax_joint.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
	
	f.ax_joint.set_xlim(1.4,2.0)
	f.ax_joint.set_ylim(-2.5,1.0)
	g.ax_joint.set_xlim(1.4,2.0)
	g.ax_joint.set_ylim(-2.5,1.0)


	g.fig.set_figwidth(8)
	g.fig.set_figheight(8)
	f.fig.set_figwidth(8)
	f.fig.set_figheight(8)

	legend_elements = [Line2D([0], [1], ls='-', color='k',label='Field'),
                       Line2D([0], [1], ls='--', color='k', label='Filament')]

	f.ax_joint.legend(handles=legend_elements,fontsize=14)

	legend_elements = [Line2D([0], [1], ls='-', color='k',label='Group'),
                       Line2D([0], [1], ls='--', color='k', label='Filament')]

	g.ax_joint.legend(handles=legend_elements,fontsize=14)

	
	g.savefig(Path+'kde_'+yprop+xprop+'_group_filament.png',eps=600)
	f.savefig(Path+'kde_'+yprop+xprop+'_field_filament.png',eps=600)
	

	plt.show()


def Binned_dist(prop,xlabel,logx):

		All_data = pd.read_csv('All_Data_Stacked.csv')
		fig = plt.figure(figsize = (6,6))
		ax = fig.gca()
   
		if logx=='y':

			All_data['log'+prop] = np.log10(All_data[prop])
		
			#Binning:
			#Binning age:
			cut_points = [0.5,1.0]

			All_data["d_Bin"] = binning(All_data["d_per"], cut_points)
			#print pd.value_counts(All_data["d_Bin"], sort=True)

			bin1 = All_data[All_data["d_Bin"]==0]
			bin1 = bin1[bin1['SFR']!=0.0]

			bin2 = All_data[All_data["d_Bin"]==1]
			bin2 = bin2[bin2['SFR']!=0.0]

			bin3 = All_data[All_data["d_Bin"]==2]
			bin3 = bin3[bin3['SFR']!=0.0]

			



			bin1['log'+prop].plot.hist(ax=ax,histtype=u'step',density=1,bins=20,label='<0.5 Mpc')
			bin2['log'+prop].plot.hist(ax=ax,histtype=u'step',density=1,bins=20,label='0.5-1.0 Mpc')
			bin3['log'+prop].plot.hist(ax=ax,histtype=u'step',density=1,bins=20,label='>1.0 Mpc')
			#bin4['log'+prop].plot.hist(ax=ax,histtype=u'step',density=1,label='0.6 - 0.8 Mpc')

			ax.set_xlabel(r'$\mathrm{'+xlabel+'}$')

		
		plt.legend(title=r'$\mathrm{d_{per}}$')
		fig.savefig(Path+'/Binned_dist_'+prop+'.png',dpi=600)

		plt.show()


def Binned_Mass():

		All_data = pd.read_csv('All_Data_Stacked.csv')
		fig = plt.figure(figsize = (6,6))
		ax = fig.gca()
   
		
		All_data['logSM'] = np.log10(All_data['M'])
		All_data['g_minus_r'] = All_data['g'] - All_data['r']

		#Binning:
		#Binning age:
		cut_points = [9.5]

		All_data = All_data[All_data['ENV'] == 'filament']

		All_data.sort_values(by='logSM',axis=0,inplace=True)

		#print All_data['logSM']

		All_data["M_Bin"] = binning(All_data["logSM"], cut_points)
		#print pd.value_counts(All_data["d_Bin"], sort=True)

		bin1 = All_data[All_data["M_Bin"]==0]
		bin2 = All_data[All_data["M_Bin"]==1]
		#bin3 = All_data[All_data["M_Bin"]==2]

		print bin1.shape,bin2.shape

		feature=['d_per','g_minus_r']

		stack_Data1 = bin1[feature]
		stack_Data1 = stack_Data1.dropna()
		stack_Data1.sort_values('d_per',inplace=True)
		stack_Data1['D'] =  stack_Data1.rolling(6000).median()['d_per']
		stack_Data1['col'] =  stack_Data1.rolling(6000).median()['g_minus_r']
		stack_Data1.plot(ax = ax,x='D',y='col',legend=False,linewidth=1.5,color='r',label=r'< 9')

		stack_Data2 = bin2[feature]
		stack_Data2 = stack_Data2.dropna()
		stack_Data2.sort_values('d_per',inplace=True)
		stack_Data2['D'] =  stack_Data2.rolling(1000).median()['d_per']
		stack_Data2['col'] =  stack_Data2.rolling(1000).median()['g_minus_r']
		stack_Data2.plot(ax = ax,x='D',y='col',legend=False,linewidth=1.5,color='b',label=r'> 9')
		
		'''
		stack_Data3 = bin3[feature]
		stack_Data3 = stack_Data3.dropna()
		stack_Data3.sort_values('d_per',inplace=True)
		stack_Data3['D'] =  stack_Data3.rolling(50).median()['d_per']
		stack_Data3['col'] =  stack_Data3.rolling(50).median()['g_minus_r']
		stack_Data3.plot(ax = ax,x='D',y='col',legend=False,linewidth=1.5,color='g',label=r'10-11')

		'''




		ax.set_xlabel(r'$\mathrm{d_{per}}$(Mpc)',fontsize=16)
		ax.set_ylabel(r'$\mathrm{g - r}$',fontsize=16)
		ax.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
		ax.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)

	
		
		
		plt.legend(title=r"$\mathrm{log \ M_{*}}$")
		#fig.savefig(Path+'/Binned_dist_'+prop+'.png',dpi=600)

		plt.show()

def Gas_plots():
	
	fig3,ax3 = plt.subplots(1,1,figsize=(6,6))
	fig4,ax4 = plt.subplots(1,1,figsize=(6,6))



	Groups,Filament,Field = Separate_ENV()


	Filament = Filament[Filament['GM']!=0.0]
	Field = Field[Field['GM']!=0.0]
	Groups = Groups[Groups['GM']!=0.0]


	Filament['logGM'] = np.log10(Filament['GM'])

	Filament['GMF'] = Filament['logGM'] - np.log10(Filament['Tot_Mass'])
	Filament['SMF'] = Filament['logSM'] - np.log10(Filament['Tot_Mass'])
	Filament['mu'] = np.log10(Filament['GM']/(Filament['GM']+Filament['SM']))
	Groups['mu'] = np.log10(Groups['GM']/(Groups['GM']+Groups['SM']))
	Field['mu'] = np.log10(Field['GM']/(Field['GM']+Field['SM']))



	feature=['d_per','GMF','SMF']

	stack_Data1 = Filament[feature]
	stack_Data1 = stack_Data1.dropna()
	stack_Data1.sort_values('d_per',inplace=True)
	stack_Data1['D'] =  stack_Data1.rolling(6000).median()['d_per']
	stack_Data1['col'] =  stack_Data1.rolling(6000).median()['GMF']
	stack_Data1.plot(ax = ax3,x='D',y='col',legend=False,linestyle='-',linewidth=1.5,color='b',label=r'$\mathrm{M_{gas}/M_{tot}}$')


	stack_Data1['D'] =  stack_Data1.rolling(6000).median()['d_per']
	stack_Data1['col2'] =  stack_Data1.rolling(6000).median()['SMF']
	stack_Data1.plot(ax = ax3,x='D',y='col2',legend=False,linestyle='-',linewidth=1.5,color='r',label=r'$\mathrm{M_{*}/M_{tot}}$')

	ax3.set_ylabel(r'$\mathrm{M/M_{tot}}$',fontsize=14)
	ax3.set_xlabel(r'$\mathrm{d_{per} \ (Mpc)}$',fontsize=14)
	ax3.tick_params(axis='both', which='minor', length=3, width=2, labelsize=16)
	ax3.tick_params(axis='both', which='major', length=5, width=2, labelsize=16)


	ax3.legend(fontsize=14)

	fig3.tight_layout()

	fig3.savefig(Path+'/d_per_vs_mass_frac.png',dpi=600)




	feature=['logSM','mu']
	stack_Data1 = Filament[feature]
	stack_Data1 = stack_Data1.dropna()
	stack_Data1.sort_values('logSM',inplace=True)
	stack_Data1['D'] =  stack_Data1.rolling(6000).median()['logSM']
	stack_Data1['col'] =  stack_Data1.rolling(6000).median()['mu']
	stack_Data1.plot(ax = ax4,x='D',y='col',legend=False,linestyle='-',linewidth=1.5,color='b',label='Filament')


	feature=['logSM','mu']
	stack_Data1 = Groups[feature]
	stack_Data1 = stack_Data1.dropna()
	stack_Data1.sort_values('logSM',inplace=True)
	stack_Data1['D'] =  stack_Data1.rolling(6000).median()['logSM']
	stack_Data1['col'] =  stack_Data1.rolling(6000).median()['mu']
	stack_Data1.plot(ax = ax4,x='D',y='col',legend=False,linestyle='-',linewidth=1.5,color='r',label='Groups')

	
	feature=['logSM','mu']
	stack_Data1 = Field[feature]
	stack_Data1 = stack_Data1.dropna()
	stack_Data1.sort_values('logSM',inplace=True)
	stack_Data1['D'] =  stack_Data1.rolling(6000).median()['logSM']
	stack_Data1['col'] =  stack_Data1.rolling(6000).median()['mu']
	stack_Data1.plot(ax = ax4,x='D',y='col',legend=False,linestyle='-',linewidth=1.5,color='g',label='Field')

	

	ax4.set_ylabel(r'$\mathrm{log_{10}(\mu)}$',fontsize=14)
	ax4.set_xlabel(r'$\mathrm{log_{10} \ M_{*}(M_{*})}$',fontsize=14)
	ax4.tick_params(axis='both', which='minor', length=3, width=2, labelsize=16)
	ax4.tick_params(axis='both', which='major', length=5, width=2, labelsize=16)


	ax4.legend(fontsize=14)

	fig4.tight_layout()


	plt.show()


def random_d_per():


		fig,ax = plt.subplots(1,1,figsize=(8,8))

		Total_data = pd.read_csv('/home/ankit/Fortran_code/All_Data_Stacked.csv')

		Filament = Total_data[Total_data['ENV']=='filament']

		d_per = Filament['d_per']

		#print Filament['d_per'].head()


		d_per = d_per.sample(frac=1)

		#print d_per.head()

		d_per.reset_index(inplace=True, drop=True)

		Filament['d_per'] = d_per

		#print Filament['d_per'].head()

		Filament.sort_values('d_per',inplace=True)

		Filament = Filament.reset_index(drop=True)

		Filament['g_minus_r'] = Filament['g'] - Filament['r']

		#print Filament['d_per'].max()


		feature=['d_per','g_minus_r']

		stack_Data = Filament[feature]


		stack_Data = stack_Data.dropna()
		stack_Data.sort_values('d_per',inplace=True)

		stack_Data.to_csv('Stacked_Data.csv',index=False)

		stack_Data['D'] =  stack_Data.rolling(6000).median()['d_per']
		stack_Data['col'] =  stack_Data.rolling(6000).median()['g_minus_r']


		stack_Data.plot(ax = ax,x='D',y='col',legend=False,linewidth=1.5,color='k')

		ax.set_xlabel(r'$\mathrm{d_{per}}$(Mpc)',fontsize=16)
		ax.set_ylabel(r'$\mathrm{g - r}$',fontsize=16)
		ax.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
		ax.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)

		#ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

		fig.tight_layout()
		fig.savefig('/home/ankit/Fortran_code/Figures/Random_dper.png',dpi=600)
		plt.show()


def Metallicity_plot():

		Total_data = pd.read_csv('/home/ankit/Fortran_code/All_Data_Stacked.csv')


		Colors = ['#3F6C45','#CB6318','#29a8ab']

		Groups,Filament,Field = Separate_ENV()


		fig,(ax1,ax2) = plt.subplots(2,1,figsize=(6,6),sharex=True)


		Filament.sort_values('SM',inplace=True)
		Filament = Filament.reset_index(drop=True)
		Filament = Filament[Filament['SF_Metal']!=0.0]
		Filament['logSM'] = np.log10(Filament['SM'])
		Filament['logSF_Metal'] = np.log10(Filament['SF_Metal'])

		feature=['logSM','logSF_Metal']
		stack_Data = Filament[feature]
		stack_Data = stack_Data.dropna()
		stack_Data.sort_values('logSM',inplace=True)
		stack_Data['D'] =  stack_Data.rolling(6000).median()['logSM']
		stack_Data['col'] =  stack_Data.rolling(6000).median()['logSF_Metal']
		stack_Data.plot(ax = ax1,x='D',y='col',legend=False,linewidth=1.5,color=Colors[1],label='Filament')

		ax1.set_xlabel(r'$\mathrm{log \ M_{*}\; (M_{\odot})}$',fontsize=16)
		ax1.set_ylabel(r'$\mathrm{Z_{SF}}$',fontsize=16)
		ax1.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
		ax1.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)



		Groups.sort_values('SM',inplace=True)
		Groups = Groups.reset_index(drop=True)
		Groups = Groups[Groups['SF_Metal']!=0.0]
		Groups['logSM'] = np.log10(Groups['SM'])
		Groups['logSF_Metal'] = np.log10(Groups['SF_Metal'])

		feature=['logSM','logSF_Metal']
		stack_Data = Groups[feature]
		stack_Data = stack_Data.dropna()
		stack_Data.sort_values('logSM',inplace=True)
		stack_Data['D'] =  stack_Data.rolling(6000).median()['logSM']
		stack_Data['col'] =  stack_Data.rolling(6000).median()['logSF_Metal']
		stack_Data.plot(ax = ax1,x='D',y='col',legend=False,linewidth=1.5,color=Colors[0],label='Group')

		ax1.set_xlabel(r'$\mathrm{log \ M_{*}\; (M_{\odot})}$',fontsize=16)
		ax1.set_ylabel(r'$\mathrm{Z_{SF}}$',fontsize=16)
		ax1.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
		ax1.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)



		Field.sort_values('SM',inplace=True)
		Field = Field.reset_index(drop=True)
		Field = Field[Field['SF_Metal']!=0.0]
		Field['logSM'] = np.log10(Field['SM'])
		Field['logSF_Metal'] = np.log10(Field['SF_Metal'])

		feature=['logSM','logSF_Metal']
		stack_Data = Field[feature]
		stack_Data = stack_Data.dropna()
		stack_Data.sort_values('logSM',inplace=True)
		stack_Data['D'] =  stack_Data.rolling(6000).median()['logSM']
		stack_Data['col'] =  stack_Data.rolling(6000).median()['logSF_Metal']
		stack_Data.plot(ax = ax1,x='D',y='col',legend=False,linewidth=1.5,color=Colors[2],label='Field')

		ax1.set_xlabel(r'$\mathrm{log \ M_{*}\; (M_{\odot})}$',fontsize=16)
		ax1.set_ylabel(r'$\mathrm{log \ Z_{SF}}$',fontsize=16)
		ax1.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
		ax1.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
		ax1.locator_params(axis='y', nbins=3)


		ax1.legend(fontsize=14)
		#ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
		#fig.tight_layout()
		#fig.savefig('/home/ankit/Fortran_code/Figures/SFMetallicity_vs_Mstar.png',dpi=600)
		#plt.show()




		#fig,ax = plt.subplots(2,1,figsize=(6,6))
		Groups,Filament,Field = Separate_ENV()

		Filament.sort_values('SM',inplace=True)
		Filament = Filament.reset_index(drop=True)
		Filament = Filament[Filament['NSF_Metal']!=0.0]
		Filament['logSM'] = np.log10(Filament['SM'])
		Filament['logNSF_Metal'] = np.log10(Filament['NSF_Metal'])

		feature=['logSM','logNSF_Metal']
		stack_Data = Filament[feature]
		stack_Data = stack_Data.dropna()
		stack_Data.sort_values('logSM',inplace=True)
		stack_Data['D'] =  stack_Data.rolling(6000).median()['logSM']
		stack_Data['col'] =  stack_Data.rolling(6000).median()['logNSF_Metal']
		stack_Data.plot(ax = ax2,x='D',y='col',legend=False,linewidth=1.5,color=Colors[1],label='Filament')

		Groups.sort_values('SM',inplace=True)
		Groups = Groups.reset_index(drop=True)
		Groups = Groups[Groups['NSF_Metal']!=0.0]
		Groups['logSM'] = np.log10(Groups['SM'])
		Groups['logNSF_Metal'] = np.log10(Groups['NSF_Metal'])

		feature=['logSM','logNSF_Metal']
		stack_Data = Groups[feature]
		stack_Data = stack_Data.dropna()
		stack_Data.sort_values('logSM',inplace=True)
		stack_Data['D'] =  stack_Data.rolling(6000).median()['logSM']
		stack_Data['col'] =  stack_Data.rolling(6000).median()['logNSF_Metal']
		stack_Data.plot(ax = ax2,x='D',y='col',legend=False,linewidth=1.5,color=Colors[0],label='Group')

		Field.sort_values('SM',inplace=True)
		Field = Field.reset_index(drop=True)
		Field = Field[Field['NSF_Metal']!=0.0]
		Field['logSM'] = np.log10(Field['SM'])
		Field['logNSF_Metal'] = np.log10(Field['NSF_Metal'])

		feature=['logSM','logNSF_Metal']
		stack_Data = Field[feature]
		stack_Data = stack_Data.dropna()
		stack_Data.sort_values('logSM',inplace=True)
		stack_Data['D'] =  stack_Data.rolling(6000).median()['logSM']
		stack_Data['col'] =  stack_Data.rolling(6000).median()['logNSF_Metal']
		stack_Data.plot(ax = ax2,x='D',y='col',legend=False,linewidth=1.5,color=Colors[2],label='Field')

		ax2.set_xlabel(r'$\mathrm{log \ M_{*}\; (M_{\odot})}$',fontsize=16)
		ax2.set_ylabel(r'$\mathrm{log \ Z_{NSF}}$',fontsize=16)
		ax2.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
		ax2.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)


		#ax2.legend(fontsize=14)
		ax2.locator_params(axis='y', nbins=3)
		
		#ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
		fig.tight_layout()
		fig.subplots_adjust(wspace=0, hspace=0)

		fig.savefig('/home/ankit/Fortran_code/Figures/Z_plot1.png',dpi=600)
		#plt.show()




		fig,(ax1,ax2) = plt.subplots(2,1,figsize=(6,6),sharex=True)

		Total_data = pd.read_csv('/home/ankit/Fortran_code/All_Data_Stacked.csv')

		Total_data.sort_values('d_per',inplace=True)
		Total_data = Total_data.reset_index(drop=True)
		Total_data = Total_data[Total_data['NSF_Metal']!=0.0]
		#Filament['logSM'] = np.log10(Filament['SM'])
		Total_data['logNSF_Metal'] = np.log10(Total_data['NSF_Metal'])

		feature=['d_per','logNSF_Metal']
		stack_Data = Total_data[feature]
		stack_Data = stack_Data.dropna()
		stack_Data.sort_values('d_per',inplace=True)
		stack_Data['D'] =  stack_Data.rolling(6000).median()['d_per']
		stack_Data['col'] =  stack_Data.rolling(6000).median()['logNSF_Metal']
		stack_Data.plot(ax = ax2,x='D',y='col',legend=False,linewidth=1.5,color='k',label='Filament')


		ax2.set_xlabel(r'$\mathrm{d_{per}}$(Mpc)',fontsize=16)
		ax2.set_ylabel(r'$\mathrm{log \ Z_{NSF}}$',fontsize=16)
		ax2.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
		ax2.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)

		#fig.savefig('/home/ankit/Fortran_code/Figures/NSFMetallicity_vs_d_per.png',dpi=600)


		#fig,ax = plt.subplots(1,1,figsize=(6,6))

		Total_data = pd.read_csv('/home/ankit/Fortran_code/All_Data_Stacked.csv')

		Total_data.sort_values('d_per',inplace=True)
		Total_data = Total_data.reset_index(drop=True)
		Total_data = Total_data[Total_data['SF_Metal']!=0.0]
		#Filament['logSM'] = np.log10(Filament['SM'])
		Total_data['logSF_Metal'] = np.log10(Total_data['SF_Metal'])

		feature=['d_per','logSF_Metal']
		stack_Data = Total_data[feature]
		stack_Data = stack_Data.dropna()
		stack_Data.sort_values('d_per',inplace=True)
		stack_Data['D'] =  stack_Data.rolling(6000).median()['d_per']
		stack_Data['col'] =  stack_Data.rolling(6000).median()['logSF_Metal']
		stack_Data.plot(ax = ax1,x='D',y='col',legend=False,linewidth=1.5,color='k',label='Filament')


		ax1.set_xlabel(r'$\mathrm{d_{per}}$(Mpc)',fontsize=16)
		ax1.set_ylabel(r'$\mathrm{log \ Z_{SF}}$',fontsize=16)
		ax1.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
		ax1.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
		ax1.locator_params(axis='y', nbins=3)
		
		fig.tight_layout()
		fig.subplots_adjust(wspace=0, hspace=0)

		fig.savefig('/home/ankit/Fortran_code/Figures/Z_plot2.png',dpi=600)
		plt.show()


def Contour_Active(props,xlabel,ylabel,logx,logy):
	
	#import matplotlib.patches as mpatches
	from matplotlib.lines import Line2D

	import seaborn as sns; sns.set(color_codes=True)
	sns.set_context(context='paper')
	sns.set_style("white")

	#plt.style.use("dark_background")

	#sns.set_style("white")

	new_plt = sns.light_palette("green",as_cmap=True)
	plt_fil = sns.dark_palette("purple", as_cmap=True)


	xprop,yprop = props

	Groups,Filament,Field = Separate_ENV()
	Groups = Groups.loc[Groups['State']=='Active']
	Filament = Filament.loc[Filament['State']=='Active']
	Field = Field.loc[Field['State']=='Active']

	Groups['log'+xprop] = np.log10(Groups[xprop])
	Filament['log'+xprop] = np.log10(Filament[xprop])
	Field['log'+xprop] = np.log10(Field[xprop])

	Groups['log'+yprop] = np.log10(Groups[yprop])
	Filament['log'+yprop] = np.log10(Filament[yprop])
	Field['log'+yprop] = np.log10(Field[yprop])


	g = sns.JointGrid(x='log'+xprop, y='log'+yprop, data=Groups,height=6)
	
	sns.kdeplot(
	    Filament['log'+xprop], Filament['log'+yprop], cmap="Greys",
	    shade=False, shade_lowest=False, ax=g.ax_joint,
	    linestyles='--',label='Field',n_levels=5
	)
	
	sns.kdeplot(
	    Groups['log'+xprop], Groups['log'+yprop], cmap="Greys",
	    shade=False, shade_lowest=False, ax=g.ax_joint,
	    linestyles='-',label='Group',n_levels=5
	)
	

	sns.distplot(
	    Groups['log'+xprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls ='-'), ax=g.ax_marg_x,axlabel=False
	)
	sns.distplot(
	    Filament['log'+xprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls='--'), ax=g.ax_marg_x,axlabel=False
	)
	sns.distplot(
	    Groups['log'+yprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls ='-'), ax=g.ax_marg_y, vertical=True,axlabel=False
	)
	sns.distplot(
	    Filament['log'+yprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls='--'), ax=g.ax_marg_y, vertical=True,axlabel=False
	)
	
	#legend_properties = {'size':12,'color':['b','r']}
	#legendMain=g.ax_joint.legend(prop=legend_properties,loc='upper right')

	#patch1 = mpatches.Patch(color='r', label='Group')
	#patch2 = mpatches.Patch(color='b', label='Filament')

	#all_handles = (patch1, patch2)

	#leg = g.ax_joint.legend(handles=all_handles)
	#g.ax_joint.add_artist(leg)

	g.ax_joint.set_xlabel(r"$\mathrm{"+xlabel+"}$",fontsize=14)
	g.ax_joint.set_ylabel(r"$\mathrm{"+ylabel+"}$",fontsize=14)

	g.ax_joint.minorticks_on()
	g.ax_joint.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	g.ax_joint.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
	
	

	######################################################################

	f = sns.JointGrid(x='log'+xprop, y='log'+yprop, data=Filament,height=6)
	
	sns.kdeplot(
	    Filament['log'+xprop], Filament['log'+yprop], cmap="Greys",
	    shade=False, shade_lowest=False, ax=f.ax_joint,
	    linestyles='--',label='Filament',n_levels=5
	)
	sns.kdeplot(
	    Field['log'+xprop], Field['log'+yprop], cmap="Greys",
	    shade=False, shade_lowest=False, ax=f.ax_joint,
	    linestyles='-',label='Field',n_levels=5
	)
	sns.distplot(
	    Filament['log'+xprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls ='--'), ax=f.ax_marg_x,axlabel=False
	)
	sns.distplot(
	    Field['log'+xprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls='-'), ax=f.ax_marg_x,axlabel=False
	)
	sns.distplot(
	    Filament['log'+yprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls ='--'), ax=f.ax_marg_y, vertical=True,axlabel=False
	)
	sns.distplot(
	    Field['log'+yprop], kde=True, hist=False, color="k",
	    kde_kws=dict(ls='-'), ax=f.ax_marg_y, vertical=True,axlabel=False
	)
	
	#legend_properties = {'size':12,'color':['b','r']}
	#legendMain=f.ax_joint.legend(prop=legend_properties,loc='upper right')

	#patch2 = mpatches.Patch(color='b', label='Filament')
	#patch1 = mpatches.Patch(color='r', label='Field')

	#all_handles = (patch1, patch2)

	#leg = f.ax_joint.legend(handles=all_handles)
	#f.ax_joint.add_artist(leg)

	f.ax_joint.set_xlabel(r"$\mathrm{"+xlabel+"}$",fontsize=14)
	f.ax_joint.set_ylabel(r"$\mathrm{"+ylabel+"}$",fontsize=14)
	
	f.ax_joint.minorticks_on()
	f.ax_joint.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	f.ax_joint.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
	
	#f.ax_joint.set_xlim(8.0,11.0)
	f.ax_joint.set_ylim(-2.5,1.0)
	#g.ax_joint.set_xlim(8.0,11.0)
	g.ax_joint.set_ylim(-2.5,1.0)


	g.fig.set_figwidth(8)
	g.fig.set_figheight(8)
	f.fig.set_figwidth(8)
	f.fig.set_figheight(8)

	legend_elements = [Line2D([0], [1], ls='-', color='k',label='Field'),
                       Line2D([0], [1], ls='--', color='k', label='Filament')]

	f.ax_joint.legend(handles=legend_elements,fontsize=14)

	legend_elements = [Line2D([0], [1], ls='-', color='k',label='Group'),
                       Line2D([0], [1], ls='--', color='k', label='Filament')]

	g.ax_joint.legend(handles=legend_elements,fontsize=14)

	
	g.savefig(Path+'Active_'+'kde_'+yprop+xprop+'_group_filament.png',dpi=1000)
	f.savefig(Path+'Active_'+'kde_'+yprop+xprop+'_field_filament.png',dpi=1000)
	

	plt.show()


def Passive_Active_Plot():

		All_data = pd.read_csv('/home/ankit/Fortran_code/All_Data_Stacked.csv')

		All_data['State'] = 'Active'

		All_data['logSM'] = np.log10(All_data['SM'])
		All_data = All_data.loc[All_data['GM']>0.0]

		All_data['logGM'] = np.log10(All_data['GM'])
		All_data['logsigma'] = np.log10(All_data['VelDisp'])

		All_data = All_data.loc[All_data['SFR']>0.0]
		All_data['logSFR'] = np.log10(All_data['SFR'])

		All_data.loc[All_data['d_per'] > 1.0,'ENV'] = 'field'

		All_data.loc[All_data['SFR']*1.0e9/All_data['SM'] <= 0.01, 'State'] = 'Passive'

		import matplotlib.pyplot as plt

		import seaborn as sns
		sns.set(style="ticks")


		#g = sns.pairplot(All_data, hue="State",vars=['logSM','logGM','logSFR'],plot_kws={'legend':False})
		g = sns.PairGrid(All_data, hue="State",vars=['logSM','logGM','logSFR'])
		g = g.map_diag(sns.kdeplot,shade=True)
		g = g.map_offdiag(plt.scatter,s=8,alpha=1.0)


		for i, j in zip(*np.triu_indices_from(g.axes, 1)):
		    g.axes[i, j].set_visible(False)

		g.axes[0,0].yaxis.set_label_text(r'$\mathrm{log \ M_{*} \ (M_{\odot})}$',fontsize=14)
		g.axes[1,0].yaxis.set_label_text(r'$\mathrm{log \ M_{gas} \ (M_{\odot})}$',fontsize=14)
		g.axes[2,0].yaxis.set_label_text(r'$\mathrm{log \ SFR \ (M_{\odot} yr^{-1})}$',fontsize=14)

		g.axes[2,0].xaxis.set_label_text(r'$\mathrm{log \ M_{*} \ (M_{\odot})}$',fontsize=14)
		g.axes[2,1].xaxis.set_label_text(r'$\mathrm{log \ M_{gas} \ (M_{\odot})}$',fontsize=14)
		g.axes[2,2].xaxis.set_label_text(r'$\mathrm{log \ SFR \ (M_{\odot} yr^{-1})}$',fontsize=14)



		labels = g._legend_data.keys()
		g.fig.legend(fontsize=14,labels=labels,title='State', loc='upper center',markerscale=4.)
		
		plt.savefig(Path+'/Passive_Active_Plot.png',dpi=600)
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
	print "\n 7. scatter plot \n"
	print "\n 8. Contour plot \n"
	print "\n 9. Mass distribution for binned d_per \n"
	print "\n 10. M hist for d binned \n"
	print "\n 11. Plots for Gas Mass \n"
	print "\n 12. Plots for Metallicity \n"
	print "\n 13. Plots for Active Contour \n"
	print "\n 14. Plots for Active Passive \n"


	opt=input('option?')

	if opt==0:
		All_data()
	if opt==1:
		Mvscol()
	if opt==2:
		prop = ConfigSectionMap("histogram")['prop']
		logit = ConfigSectionMap("histogram")['log']
		dens = ConfigSectionMap("histogram")['dens']
		histogram(prop,logit,dens)

	if opt==3:
		prop1 = ConfigSectionMap("JointPlot")['prop1']

		prop2 = ConfigSectionMap("JointPlot")['prop2']

		props = (prop2,prop1)
		
		xlabel=ConfigSectionMap("JointPlot")['xlabel']
		
		ylabel=ConfigSectionMap("JointPlot")['ylabel']

		remove0x=ConfigSectionMap("JointPlot")['remove0x']

		remove0y=ConfigSectionMap("JointPlot")['remove0y']

		logx=ConfigSectionMap("JointPlot")['logx']

		logy=ConfigSectionMap("JointPlot")['logy']

		
		jointplot(props,xlabel,ylabel,remove0x,remove0y,logx,logy)

	if opt==4:
		Voilin_plot()
	if opt==5:

		sliceno=ConfigSectionMap("Colslice")['sliceno']

		Colored_slice(sliceno)

	if opt==6:
		dvscol()

	if opt==7:

		prop1 = ConfigSectionMap("Scatter")['prop1']

		prop2 = ConfigSectionMap("Scatter")['prop2']

		props = (prop2,prop1)

		colprop=ConfigSectionMap("Scatter")['colprop']

		xlabel=ConfigSectionMap("Scatter")['xlabel']
		
		ylabel=ConfigSectionMap("Scatter")['ylabel']

		logx=ConfigSectionMap("Scatter")['logx']
		logy=ConfigSectionMap("Scatter")['logy']
		logcon=ConfigSectionMap("Scatter")['logcol']

		Scatter(props,colprop,xlabel,ylabel,logx,logy,logcon)

	if opt==8:

		prop1 = ConfigSectionMap("Contour")['prop1']
		prop2 = ConfigSectionMap("Contour")['prop2']
		props = (prop2,prop1)
		xlabel=ConfigSectionMap("Contour")['xlabel']
		ylabel=ConfigSectionMap("Contour")['ylabel']
		logx=ConfigSectionMap("Contour")['logx']
		logy=ConfigSectionMap("Contour")['logy']
		
		Contour(props,xlabel,ylabel,logx,logy)
	
	if opt==9:

		prop = ConfigSectionMap("Binned")['prop']
		xlabel=ConfigSectionMap("Binned")['xlabel']
		logx=ConfigSectionMap("Binned")['logx']

		Binned_dist(prop,xlabel,logx)

	if opt==10:

		Binned_Mass()
	
	if opt==11:

		Gas_plots()

	if opt==12:

		Metallicity_plot()

	if opt==13:

		prop1 = ConfigSectionMap("Contour_Active")['prop1']
		prop2 = ConfigSectionMap("Contour_Active")['prop2']
		props = (prop2,prop1)
		xlabel=ConfigSectionMap("Contour_Active")['xlabel']
		ylabel=ConfigSectionMap("Contour_Active")['ylabel']
		logx=ConfigSectionMap("Contour_Active")['logx']
		logy=ConfigSectionMap("Contour_Active")['logy']
		
		Contour_Active(props,xlabel,ylabel,logx,logy)
	
	if opt==14:

		Passive_Active_Plot()

if __name__=='__main__':

	main()
