import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pylab import savefig
import numpy as np

import pickle
import warnings
from matplotlib.pyplot import figure

title=f'Selection Pressures of 10,002nt sequence p={0.001}'

wd='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/fmh_dnds_not_using_sourmash_translate'

positive_fmh_dnds_folder1="/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/fmh_dnds/positive_selection"
negative_fmh_dnds_folder1="/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/fmh_dnds/negative_selection"
fmhdnds1='all_dnds_constant.csv'

positive_fmh_dnds_folder="/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/fmh_dnds_not_using_sourmash_translate/positive_selection"
negative_fmh_dnds_folder="/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/fmh_dnds_not_using_sourmash_translate/negative_selection"
fmhdnds='dnds_constant_all.csv'

kaks_folder="/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/kaks_NG"
negative_kaks='negative_selection_queries_10002_0.001.axt.kaks'
positive_kaks='positive_selection_queries_10002_0.001.axt.kaks'

#read in pivot tables
methods = ['NG']
kaks_calc = pd.read_csv(f'{kaks_folder}/{negative_kaks}',sep='\t')
kaks_calc = kaks_calc.pivot(index='Sequence',columns='Method')[['Ka/Ks']] #use when no duplicates
kaks_calc = kaks_calc['Ka/Ks'][methods].reset_index().rename(columns={'NG': 'negative-NG'})
kaks_calc['Sequence'] = kaks_calc['Sequence'].str.replace(r'negative_', r'gene_')
print(kaks_calc)

methods = ['NG']
kaks_calc2 = pd.read_csv(f'{kaks_folder}/{positive_kaks}',sep='\t')
kaks_calc2 = kaks_calc2.pivot(index='Sequence',columns='Method')[['Ka/Ks']] #use when no duplicates
kaks_calc2 = kaks_calc2['Ka/Ks'][methods].reset_index().rename(columns={'NG': 'positive-NG'})
kaks_calc2['Sequence'] = kaks_calc2['Sequence'].str.replace(r'positive_', r'gene_')
print(kaks_calc2)

#file is contanation of all dnds_contant files but make sure headers are only in the first line of the file
ksizes = ['5','7']
fmh_dnds = pd.read_csv(f'{positive_fmh_dnds_folder}/{fmhdnds}',sep=',')
fmh_dnds = fmh_dnds[(fmh_dnds['A'] == 'ref_gene')].rename(columns={'ksize': 'Method','B':'Sequence'}).pivot(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
fmh_dnds = fmh_dnds['dNdS_ratio_constant'][ksizes].rename(columns={'5': 'positive-5','7':'positive-7'}).reset_index()
fmh_dnds['Sequence'] = fmh_dnds['Sequence'].str.replace(r'positive_', r'gene_')


#file is contanation of all dnds_contant files but make sure headers are only in the first line of the file
fmh_dnds2 = pd.read_csv(f'{negative_fmh_dnds_folder}/{fmhdnds}',sep=',')
fmh_dnds2 = fmh_dnds2[(fmh_dnds2['A'] == 'ref_gene')].rename(columns={'ksize': 'Method','B':'Sequence'}).pivot(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
fmh_dnds2 = fmh_dnds2['dNdS_ratio_constant'][ksizes].rename(columns={'5': 'negative-5','7':'negative-7'}).reset_index()
fmh_dnds2['Sequence'] = fmh_dnds2['Sequence'].str.replace(r'negative_', r'gene_')

#file is contanation of all dnds_contant files but make sure headers are only in the first line of the file
ksizes = [5,7]
fmh_dnds3 = pd.read_csv(f'{positive_fmh_dnds_folder1}/{fmhdnds1}',sep=',')
fmh_dnds3 = fmh_dnds3[(fmh_dnds3['A'] == 'ref_gene')].rename(columns={'ksize': 'Method','B':'Sequence'}).pivot(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
fmh_dnds3 = fmh_dnds3['dNdS_ratio_constant'][ksizes].rename(columns={5: 'sm_translate_positive-5',7:'sm_translate_positive-7'}).reset_index()
fmh_dnds3['Sequence'] = fmh_dnds3['Sequence'].str.replace(r'positive_', r'gene_')


#file is contanation of all dnds_contant files but make sure headers are only in the first line of the file
fmh_dnds4 = pd.read_csv(f'{negative_fmh_dnds_folder1}/{fmhdnds1}',sep=',')
fmh_dnds4 = fmh_dnds4[(fmh_dnds4['A'] == 'ref_gene')].rename(columns={'ksize': 'Method','B':'Sequence'}).pivot(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
fmh_dnds4 = fmh_dnds4['dNdS_ratio_constant'][ksizes].rename(columns={5: 'sm_translate_negative-5',7:'sm_translate_negative-7'}).reset_index()
fmh_dnds4['Sequence'] = fmh_dnds4['Sequence'].str.replace(r'negative_', r'gene_')

#merge dataframes
df1 = pd.merge(kaks_calc, kaks_calc2, 'left', on = ["Sequence"])
df1 = pd.merge(df1, fmh_dnds, 'left', on = ["Sequence"])
df1 = pd.merge(df1, fmh_dnds2, 'left', on = ["Sequence"])
df1 = pd.merge(df1, fmh_dnds3, 'left', on = ["Sequence"])
df1 = pd.merge(df1, fmh_dnds4, 'left', on = ["Sequence"])
print(df1)
#methods_new=['negative-NG','sm_translate_negative-5','negative-5','sm_translate_negative-7','negative-7','positive-NG','sm_translate_positive-5','positive-5','sm_translate_positive-7','positive-7']
methods_new=['negative-NG','sm_translate_negative-5','negative-5','sm_translate_negative-7','negative-7']
df1=df1[methods_new]
#print(df1)

#filter outliers
df1 = df1.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio
df1 = df1[(df1 != 0).all(1)]
df1 = df1.sort_values(by='negative-NG')
df1 = df1.astype(float)
#print(df1.round(6))

####################### heatmap #######################

plt = sns.heatmap(df1, vmin=0, vmax=2,xticklabels=True, cmap='seismic')
plt.set_title(title)
plt.figure.savefig(f'{wd}/redo_0.001_negative_all_methods_heatmap.png',bbox_inches='tight')




