import pandas as pd
import sys
import matplotlib.pyplot as plt


numslices=int(sys.argv[1])


fig,ax = plt.subplots(1,1,figsize=(8,8))


stack_Data = pd.DataFrame()

for i in range(1,numslices+1):

	direc = './SLICED_'+ str(i)

	Total_data = pd.read_csv(direc+'/Final_data.csv')

	Total_data = Total_data[Total_data['ENV']=='filament']

	Total_data.sort_values('d_per',inplace=True)

	Total_data = Total_data.reset_index(drop=True)

	Total_data['g_minus_r'] = Total_data['g'] - Total_data['r']


	
	feature=['d_per','g_minus_r']

	sd = Total_data[feature]

	stack_Data=pd.concat([stack_Data,sd],ignore_index=True)

	

	#Mean_col = Mean_col + Total_data['D'].dropna()


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

ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

fig.savefig('d_per_vs_gminusr.png',dpi=600)
plt.show()
