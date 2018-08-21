import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np
from os import listdir
from astropy.table import Table, vstack
from astropy.io import ascii
import ConfigParser
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata

plt.style.use('dark_background')

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

			Total_data = pd.read_csv('All_Data_Stacked.csv')
			#Total_data = pd.read_csv('SLICED_13/Final_data.csv')

			
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

			fig.savefig(Path+'slice13_Stellar_mass_vs_uminusr.png',dpi=600,bbox_inches='tight')

			plt.show()

def histogram(prop,logit,dens):
			######################################################HISTOGRAM PLOT#############################################

			fig2,ax2 = plt.subplots(1,1,figsize=(8,8))

			Groups,Filament,Field = Separate_ENV()
			
			if logit=='yes':
				Groups['log'+prop] = np.log10(Groups[prop])
				Filament['log'+prop] = np.log10(Filament[prop])
				Field['log'+prop] = np.log10(Field[prop])
				prop = 'log'+prop

			if dens=='yes':
				dens1 = True
			else:
				dens1 = False


			y1,binEdges1=np.histogram(Groups[prop],bins='auto',density=dens1)
			y2,binEdges2=np.histogram(Filament[prop],bins='auto',density=dens1)
			y3,binEdges3=np.histogram(Field[prop],bins='auto',density=dens1)

			bincenters1 = 0.5*(binEdges1[1:]+binEdges1[:-1])
			bincenters2 = 0.5*(binEdges2[1:]+binEdges2[:-1])
			bincenters3 = 0.5*(binEdges3[1:]+binEdges3[:-1])

			ax2.plot(binEdges1[:-1], y1,'-o', color='m', linewidth=2, ms=8, label='Groups')
			ax2.plot(binEdges2[:-1], y2,'-o', color='peru', linewidth=2, ms=8, label='Filament')
			ax2.plot(binEdges3[:-1], y3,'-o', color='r', linewidth=2, ms=8   ,label='Field')  

			ax2.set_xlabel(prop)
			ax2.set_ylabel('Probability density')

			fig2.legend()
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
			plt.style.use('dark_background')

			fig,ax = plt.subplots(1,1,figsize=(8,8))

			Groups,Filament,Field = Separate_ENV()
			ax = sns.violinplot(x="ENV", y="u_minus_r",
			        data=Total_data, palette="muted", split=True)

			ax.set(ylabel=r'$\mathrm{u - r}$',xlabel='Environment')
			fig.savefig(Path+'violinplot.png',dpi=600)
			plt.show()


def Colored_slice():

			################################################################################################################



			FIG1,AX1 = plt.subplots(1,1,figsize=(8,8))
			global direc
			direc = './SLICED_' + str(15) 
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
			
			AX1.scatter(Field['xslice'],Field['yslice'],s=1,color=Colors[0],label='field',alpha=0.5)
			AX1.scatter(Filament['xslice'],Filament['yslice'],s=1,color=Colors[1],label='filament',alpha=0.5)
			AX1.scatter(Groups['xslice'],Groups['yslice'],s=1,color=Colors[2],label='group',alpha=0.5)

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
			
			AX1.legend(fontsize=14,markerscale=5)

		
			FIG1.savefig(Path+'/Color_without_filament.png',dpi=600)
			
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
			'Tot_Mass','SM_sh','Metal','SF_Metal','NSF_Metal','SF_O','SF_H','x','y','z','Vel','Mass_sh','u','g','r','i',\
			'zmag','Y','J','H','K','index','ENV','d_long','d_per','d_tot','xslice','yslice'])
		
		D_final = vstack([D_final, Di])

		print len(D_final['d_per'])


	#ascii.write('All_Data_Stacked.csv',D_final)
	Write_DF = D_final.to_pandas()

	Write_DF.to_csv('All_Data_Stacked.csv',index=False)


def dvscol():

	fig,ax = plt.subplots(1,1,figsize=(8,8))

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
	stack_Data.sort_values('d_per',inplace=True)

	stack_Data.to_csv('Stacked_Data.csv',index=False)

	stack_Data['D'] =  stack_Data.rolling(6000).median()['d_per']
	stack_Data['col'] =  stack_Data.rolling(6000).median()['g_minus_r']


	stack_Data.plot(ax = ax,x='D',y='col',legend=False,linewidth=1.5,color='white')

	ax.set_xlabel(r'$\mathrm{d_{per}}$(Mpc)',fontsize=16)
	ax.set_ylabel(r'$\mathrm{g - r}$',fontsize=16)
	ax.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	ax.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)

	#ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

	fig.savefig(Path+'d_per_vs_gminusr.png',dpi=600)
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
	
	
	import matplotlib.patches as mpatches

	import seaborn as sns; sns.set(color_codes=True)
	sns.set_context(context='paper')
	plt.style.use("dark_background")

	#sns.set_style("white")

	new_plt = sns.light_palette("green",as_cmap=True)
	plt_fil = sns.dark_palette("purple", as_cmap=True)


	xprop,yprop = props

	Groups,Filament,Field = Separate_ENV()

	Groups,Filament,Field = Preprocess(Groups,Filament,Field,xprop,yprop,'y','y',logx,logy)


	g = sns.JointGrid(x='log'+xprop, y='log'+yprop, data=Groups,height=6)
	
	sns.kdeplot(
	    Field['log'+xprop], Field['log'+yprop], cmap="Blues",
	    shade=False, shade_lowest=False, ax=g.ax_joint,
	    linestyles='-',label='Field'
	)
	
	sns.kdeplot(
	    Groups['log'+xprop], Groups['log'+yprop], cmap="Reds",
	    shade=False, shade_lowest=False, ax=g.ax_joint,
	    linestyles='-',label='Group'
	)
	

	sns.distplot(
	    Groups['log'+xprop], kde=True, hist=False, color="r",
	    kde_kws=dict(ls ='-'), ax=g.ax_marg_x,axlabel=False
	)
	sns.distplot(
	    Field['log'+xprop], kde=True, hist=False, color="b",
	    kde_kws=dict(ls='-'), ax=g.ax_marg_x,axlabel=False
	)
	sns.distplot(
	    Groups['log'+yprop], kde=True, hist=False, color="r",
	    kde_kws=dict(ls ='-'), ax=g.ax_marg_y, vertical=True,axlabel=False
	)
	sns.distplot(
	    Field['log'+yprop], kde=True, hist=False, color="b",
	    kde_kws=dict(ls='-'), ax=g.ax_marg_y, vertical=True,axlabel=False
	)
	
	#legend_properties = {'size':12,'color':['b','r']}
	#legendMain=g.ax_joint.legend(prop=legend_properties,loc='upper right')

	patch1 = mpatches.Patch(color='r', label='Group')
	patch2 = mpatches.Patch(color='b', label='Field')

	all_handles = (patch1, patch2)

	leg = g.ax_joint.legend(handles=all_handles)
	g.ax_joint.add_artist(leg)

	g.ax_joint.set_xlabel(r"$\mathrm{"+xlabel+"}$",fontsize=12)
	g.ax_joint.set_ylabel(r"$\mathrm{"+ylabel+"}$",fontsize=12)

	g.ax_joint.minorticks_on()
	g.ax_joint.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	g.ax_joint.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
	


	######################################################################

	f = sns.JointGrid(x='log'+xprop, y='log'+yprop, data=Filament,height=6)
	
	sns.kdeplot(
	    Filament['log'+xprop], Filament['log'+yprop], cmap="Reds",
	    shade=False, shade_lowest=False, ax=f.ax_joint,
	    linestyles='-',label='Filament'
	)
	sns.kdeplot(
	    Field['log'+xprop], Field['log'+yprop], cmap="Blues",
	    shade=False, shade_lowest=False, ax=f.ax_joint,
	    linestyles='-',label='Field'
	)
	sns.distplot(
	    Filament['log'+xprop], kde=True, hist=False, color="r",
	    kde_kws=dict(ls ='-'), ax=f.ax_marg_x,axlabel=False
	)
	sns.distplot(
	    Field['log'+xprop], kde=True, hist=False, color="b",
	    kde_kws=dict(ls='-'), ax=f.ax_marg_x,axlabel=False
	)
	sns.distplot(
	    Filament['log'+yprop], kde=True, hist=False, color="r",
	    kde_kws=dict(ls ='-'), ax=f.ax_marg_y, vertical=True,axlabel=False
	)
	sns.distplot(
	    Field['log'+yprop], kde=True, hist=False, color="b",
	    kde_kws=dict(ls='-'), ax=f.ax_marg_y, vertical=True,axlabel=False
	)
	
	#legend_properties = {'size':12,'color':['b','r']}
	#legendMain=f.ax_joint.legend(prop=legend_properties,loc='upper right')

	patch2 = mpatches.Patch(color='r', label='Filament')
	patch1 = mpatches.Patch(color='b', label='Field')

	all_handles = (patch1, patch2)

	leg = f.ax_joint.legend(handles=all_handles)
	f.ax_joint.add_artist(leg)

	f.ax_joint.set_xlabel(r"$\mathrm{"+xlabel+"}$",fontsize=12)
	f.ax_joint.set_ylabel(r"$\mathrm{"+ylabel+"}$",fontsize=12)
	
	f.ax_joint.minorticks_on()
	f.ax_joint.tick_params(axis='both', which='minor', length=3, width=2, labelsize=14)
	f.ax_joint.tick_params(axis='both', which='major', length=5, width=2, labelsize=14)
	
	g.fig.set_figwidth(8)
	g.fig.set_figheight(8)
	f.fig.set_figwidth(8)
	f.fig.set_figheight(8)


	g.savefig(Path+'kde_veldisp_group_field.png',eps=600)
	f.savefig(Path+'kde_veldisp_filament_field.png',eps=600)
	
	plt.show()


def Binned_dist(prop,xlabel,logx):

		All_data = pd.read_csv('All_Data_Stacked.csv')
		fig = plt.figure(figsize = (6,6))
		ax = fig.gca()
   
		if logx=='y':

			All_data['log'+prop] = np.log10(All_data[prop])
		
			#Binning:
			#Binning age:
			cut_points = [0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8]

			All_data["d_Bin"] = binning(All_data["d_per"], cut_points)
			#print pd.value_counts(All_data["d_Bin"], sort=True)

			bin1 = All_data[All_data["d_Bin"]==0]
			bin2 = All_data[All_data["d_Bin"]==1]
			bin3 = All_data[All_data["d_Bin"]==2]
			bin4 = All_data[All_data["d_Bin"]==3]


			bin1['log'+prop].plot.hist(ax=ax,histtype=u'step',density=1,label='<0.2 Mpc')
			bin2['log'+prop].plot.hist(ax=ax,histtype=u'step',density=1,label='0.2 - 0.4 Mpc')
			bin3['log'+prop].plot.hist(ax=ax,histtype=u'step',density=1,label='0.4 - 0.6 Mpc')
			bin4['log'+prop].plot.hist(ax=ax,histtype=u'step',density=1,label='0.6 - 0.8 Mpc')

			ax.set_xlabel(r'$\mathrm{'+xlabel+'}$')

		
		plt.legend()
		fig.savefig(Path+'/Binned_dist.png',dpi=600)

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

	opt=input('option?')

	if opt==0:
		All_data()
	if opt==1:
		Mvscol()
	if opt==2:
		prop = raw_input('what property?')
		logit = raw_input('want to take log? yes or no? \t')
		dens = raw_input('want the density? yes or no? \t')
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
		Colored_slice()
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

if __name__=='__main__':

	main()