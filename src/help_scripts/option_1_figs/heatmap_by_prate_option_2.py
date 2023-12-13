import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pylab import savefig
import numpy as np

import pickle
import warnings
from matplotlib.pyplot import figure

prate='0.001'
title=f'Compare Selection Pressures {prate} for Option 2'
folder1=f"/data/jzr5814/sourmash_dnds_estimation/tests/test/testing_option_2/{prate}"
folder3=f"/data/jzr5814/sourmash_dnds_estimation/tests/test/testing_option_2"
kaks='kaks_sequences.axt.kaks'
fmhdnds='dnds_constant_all.csv'

#read in pivot tables
#methods = ['GNG','GY-HKY','LPB','LWL','NG','YN']
methods = ['NG']

kaks_calc = pd.read_csv(f'{folder1}/{kaks}',sep='\t')
#kaks_calc = kaks_calc.pivot_table(index='Sequence',columns='Method')[['Ka/Ks']]
kaks_calc = kaks_calc.pivot(index='Sequence',columns='Method')[['Ka/Ks']] #use when no duplicates
kaks_calc = kaks_calc['Ka/Ks'][methods].reset_index().rename(columns={'NG': f'{prate}-NG'})
kaks_calc['Sequence'] = kaks_calc['Sequence'].str.replace(r'gene_0\.001_(\d+)', r'gene_\1')
#print(kaks_calc)

#file is contanation of all dnds_contant files but make sure headers are only in the first line of the file
#ksizes = [5,7,10,15,20]
ksizes = [7,10,15]
fmh_dnds = pd.read_csv(f'{folder1}/{fmhdnds}',sep=',').rename(columns={'ksize': 'Method','B':'Sequence'}).pivot(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
fmh_dnds = fmh_dnds['dNdS_ratio_constant'][ksizes].rename(columns={7: f'{prate}-7',10:f'{prate}-10',15:f'{prate}-15'}).reset_index()
fmh_dnds['Sequence'] = fmh_dnds['Sequence'].str.replace(r'gene_0\.001_(\d+)', r'gene_\1')
#print(fmh_dnds)

#merge dataframes
df1 = pd.merge(kaks_calc, fmh_dnds, 'left', on = ["Sequence"])
#df1.to_csv(f'{folder2}/dNdS_methods_pairwise_without_mean_without_evola_maxC.csv')
methods_new=[f'{prate}-NG',f'{prate}-7',f'{prate}-10',f'{prate}-15']
df1=df1[methods_new]
#print(df1)

#filter outliers
df1 = df1.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio
#print(df1)

df1 = df1[(df1 != 0).all(1)]

df1 = df1.sort_values(by=f'{prate}-NG')

print(df1)

####################### heatmap #######################

plt = sns.heatmap(df1, vmin=0, vmax=2,xticklabels=True, cmap='seismic')
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
plt.set_title(title)
#plt.figure.savefig(f'{folder3}/heatmap_{prate}_only_option_2.png',bbox_inches='tight')


