import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pylab import savefig
import numpy as np


wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_nt_genome_redo'

title='10 E coli. Strain nuc Genomes using Max Cfrac(A,B)'


#file is contanation of all dnds_contant files but make sure headers are only in the first line of the file
ksizes = [5,7,10,15,20]
fmh_dnds = pd.read_csv(f'{wd}/all_dnds_max_containments_with_selection_pressure_fixed_columns_names.csv',sep=',')
#replace strings
fmh_dnds['A'] = fmh_dnds['A'].str.replace('.name_change.fna', '')
fmh_dnds['B'] = fmh_dnds['B'].str.replace('.name_change.fna', '')
fmh_dnds['sequence_comparison'] = fmh_dnds['sequence_comparison'].str.replace('_vs_', ',').str.replace('.name_change.fna', '')
#prepare table for plotting
print(fmh_dnds['sequence_comparison'])

fmh_dnds = fmh_dnds.rename(columns={'ksize': 'Method', 'sequence_comparison': 'Sequence'}).pivot(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
print(type(fmh_dnds))

df1 = fmh_dnds['dNdS_ratio_constant'][ksizes]

#filter outliers
#print(type(df1))
df1 = df1.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio
#print(type(df1))
df1 = df1[(df1 != 0).all(1)]

#df1.to_csv(f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/{folder}/heatmap_pairwise_no_mean_no_evoladb.csv')
df1.to_csv(f'{wd}/heatmap_pairwise_all_dnds_max_containments_fixed_columns.csv')
#df1=df1[methods+ksizes+evola_dnds_method]
print(len(df1))
#df1 = df1.set_index('Sequence')

#plt = sns.heatmap(df1, cmap="crest", vmin=0, vmax=2)
#plt = sns.heatmap(df1, vmin=0, vmax=2, annot=False, linewidths=0.05, linecolor='white',yticklabels=True,xticklabels=True,cmap='tab10')
#cmap=sns.diverging_palette(20, 220, as_cmap=True)
#plt = sns.heatmap(df1, vmin=0, vmax=2,xticklabels=True, cmap='vlag')
plt = sns.heatmap(df1, vmin=0, vmax=2,xticklabels=True, cmap='seismic')
plt.set_title(title)
#plt.figure.savefig(f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/{folder}/heatmap_pairwise_mean_coolwarm.png',bbox_inches='tight')
plt.figure.savefig(f'{wd}/heatmap_pairwise_all_dnds_max_containments_fixed_columns.png',bbox_inches='tight')
