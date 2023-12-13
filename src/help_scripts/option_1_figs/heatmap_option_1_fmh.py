import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pylab import savefig
import numpy as np

import pickle
import warnings
from matplotlib.pyplot import figure

title='Compare Selection Pressures for Option 2'
folder1="/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/marine_bacteria"
fmhdnds='all_dnds_constant.csv'

"""
#read in pivot tables
methods = ['NG']

kaks_calc = pd.read_csv(f'{folder1}/{kaks}',sep='\t')
kaks_calc = kaks_calc.pivot(index='Sequence',columns='Method')[['Ka/Ks']] #use when no duplicates
kaks_calc = kaks_calc['Ka/Ks'][methods].reset_index().rename(columns={'NG': '0.001-NG'})
kaks_calc['Sequence'] = kaks_calc['Sequence'].str.replace(r'gene_0\.001_(\d+)', r'gene_\1')
#print(kaks_calc)
"""

#file is contanation of all dnds_contant files but make sure headers are only in the first line of the file
ksizes = [5,7,10,15,20]
#fmh_dnds = pd.read_csv(f'{folder1}/{fmhdnds}',sep=',').rename(columns={'ksize': 'Method','sequence_comparison':'Sequence'}).pivot(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
fmh_dnds = pd.read_csv(f'{folder1}/{fmhdnds}',sep=',').rename(columns={'ksize': 'Method','sequence_comparison':'Sequence'}).pivot(index=['A','B'],columns='Method')[['dNdS_ratio_constant']]
fmh_dnds = fmh_dnds['dNdS_ratio_constant'][ksizes]
#print(fmh_dnds)

#filter outliers
df1 = fmh_dnds.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio
df1 = df1[(df1 != 0).all(1)]
df1 = df1.sort_values(by=20)
print(df1)

####################### heatmap #######################
print(df1)
plt = sns.heatmap(df1, vmin=0, vmax=2,xticklabels=True, cmap='seismic')
plt.set_title(title)
plt.figure.savefig(f'{folder1}/heatmap_option_1_fmh.png',bbox_inches='tight')




