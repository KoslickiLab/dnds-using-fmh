import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pylab import savefig
import numpy as np


wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_genome_sketches'

#file is contanation of all dnds_contant files but make sure headers are only in the first line of the file
ksizes = [5,7,10,15,20]
fmh_dnds = pd.read_csv(f'{wd}/all_dnds_constant.csv',sep=',').rename(columns={'ksize': 'Method', 'sequence_comparison': 'Sequence'}).pivot(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
df1 = fmh_dnds['dNdS_ratio_constant'][ksizes]

#filter outliers
df1 = df1.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio
df1 = df1[(df1 != 0).all(1)]

#df1.to_csv(f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/{folder}/heatmap_pairwise_no_mean_no_evoladb.csv')
df1.to_csv(f'{wd}/heatmap_pairwise.csv')
#df1=df1[methods+ksizes+evola_dnds_method]
print(len(df1))
#df1 = df1.set_index('Sequence')

#plt = sns.heatmap(df1, cmap="crest", vmin=0, vmax=2)
#plt = sns.heatmap(df1, vmin=0, vmax=2, annot=False, linewidths=0.05, linecolor='white',yticklabels=True,xticklabels=True,cmap='tab10')
#cmap=sns.diverging_palette(20, 220, as_cmap=True)
#plt = sns.heatmap(df1, vmin=0, vmax=2,xticklabels=True, cmap='vlag')
plt = sns.heatmap(df1, vmin=0, vmax=2,xticklabels=True, cmap='seismic')
plt.set_title('CDS Genome Pairwise')
#plt.figure.savefig(f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/{folder}/heatmap_pairwise_mean_coolwarm.png',bbox_inches='tight')
plt.figure.savefig(f'{wd}/heatmap_pairwise.png',bbox_inches='tight')
