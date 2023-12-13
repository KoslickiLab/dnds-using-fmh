import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pylab import savefig
import numpy as np

import pickle
import warnings
from matplotlib.pyplot import figure

title='Compare Selection Pressures for Option 2'
folder1="/data/jzr5814/sourmash_dnds_estimation/tests/test/testing_option_2/0.001"
folder2="/data/jzr5814/sourmash_dnds_estimation/tests/test/testing_option_2/0.1"
folder3="/data/jzr5814/sourmash_dnds_estimation/tests/test/testing_option_2"
kaks='kaks_sequences.axt.kaks'
fmhdnds='dnds_constant_all.csv'

#read in pivot tables
methods = ['NG']

kaks_calc = pd.read_csv(f'{folder1}/{kaks}',sep='\t')
kaks_calc = kaks_calc.pivot(index='Sequence',columns='Method')[['Ka/Ks']] #use when no duplicates
kaks_calc = kaks_calc['Ka/Ks'][methods].reset_index().rename(columns={'NG': '0.001-NG'})
kaks_calc['Sequence'] = kaks_calc['Sequence'].str.replace(r'gene_0\.001_(\d+)', r'gene_\1')
#print(kaks_calc)

#file is contanation of all dnds_contant files but make sure headers are only in the first line of the file
ksizes = [7,10,15]
fmh_dnds = pd.read_csv(f'{folder1}/{fmhdnds}',sep=',').rename(columns={'ksize': 'Method','B':'Sequence'}).pivot(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
fmh_dnds = fmh_dnds['dNdS_ratio_constant'][ksizes].rename(columns={7: '0.001-7',10:'0.001-10',15:'0.001-15'}).reset_index()
fmh_dnds['Sequence'] = fmh_dnds['Sequence'].str.replace(r'gene_0\.001_(\d+)', r'gene_\1')
#print(fmh_dnds)

kaks_calc2 = pd.read_csv(f'{folder2}/{kaks}',sep='\t')
kaks_calc2 = kaks_calc2.pivot(index='Sequence',columns='Method')[['Ka/Ks']] #use when no duplicates
kaks_calc2 = kaks_calc2['Ka/Ks'][methods].reset_index().rename(columns={'NG': '0.1-NG'})
kaks_calc2['Sequence'] = kaks_calc2['Sequence'].str.replace(r'gene_0\.1_(\d+)', r'gene_\1')
#print(kaks_calc2)

#file is contanation of all dnds_contant files but make sure headers are only in the first line of the file
ksizes = [7,10,15]
fmh_dnds2 = pd.read_csv(f'{folder2}/{fmhdnds}',sep=',').rename(columns={'ksize': 'Method','B':'Sequence'}).pivot(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
fmh_dnds2 = fmh_dnds2['dNdS_ratio_constant'][ksizes].rename(columns={7: '0.1-7',10:'0.1-10',15:'0.1-15'}).reset_index()
fmh_dnds2['Sequence'] = fmh_dnds2['Sequence'].str.replace(r'gene_0\.1_(\d+)', r'gene_\1')
#print(fmh_dnds2)

#merge dataframes
df1 = pd.merge(kaks_calc, fmh_dnds, 'left', on = ["Sequence"])
df1 = pd.merge(df1, kaks_calc2, 'left', on = ["Sequence"])
df1 = pd.merge(df1, fmh_dnds2, 'left', on = ["Sequence"])
methods_new=['0.001-NG','0.1-NG','0.001-7','0.1-7','0.001-10','0.1-10','0.001-15','0.1-15']
df1=df1[methods_new]
#print(df1)

#filter outliers
df1 = df1.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio
df1 = df1[(df1 != 0).all(1)]
df1 = df1.sort_values(by='0.1-NG')
#print(df1)

####################### heatmap #######################

plt = sns.heatmap(df1, vmin=0, vmax=2,xticklabels=True, cmap='seismic')
plt.set_title(title)
plt.figure.savefig(f'{folder3}/heatmap_0.1_option_2.png',bbox_inches='tight')




