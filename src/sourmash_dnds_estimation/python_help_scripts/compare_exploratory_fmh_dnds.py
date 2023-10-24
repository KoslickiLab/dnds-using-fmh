import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import warnings
from matplotlib.pyplot import figure


wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_genome_sketches'

#file is contanation of all dnds_contant files but make sure headers are only in the first line of the file
ksizes = [5,7,10,15,20]
fmh_dnds = pd.read_csv(f'{wd}/all_dnds_constant.csv',sep=',')
print(fmh_dnds)

filter1 = fmh_dnds.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio
print(filter1)
filter2 = filter1[(filter1 != 0).all(1)]
print(filter2)
filter3 = filter2[(filter2["dNdS_ratio_constant"] >= 0) & (filter2["dNdS_ratio_constant"] <= 2)]

print(filter3)

df1 = filter3.rename(columns={'ksize': 'Method', 'sequence_comparison': 'Sequence'}).pivot(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
df1 = df1['dNdS_ratio_constant'][ksizes]

#merge dataframes
df1.to_csv(f'{wd}/dNdS_fmh_methods_pairwise.csv')


print(df1)

#df1=df1[ksizes]
#boxplot
bx_plotvals1, vals1, names1, xs1 = [],[],[],[]
for i, col in enumerate(df1.columns):
    vals1.append(df1[col].values)
    bx_plotvals1.append(df1[col][~np.isnan(df1[col].values)])
    names1.append(str(col)+"-mer")
    xs1.append(np.random.normal(i + 1, 0.04, df1[col].values.shape[0]))  # adds jitter to the data points - can be adjusted

fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(20, 10))

##### Set style options here #####
sns.set_style("whitegrid")  # "white","dark","darkgrid","ticks"
boxprops = dict(linestyle='-', linewidth=1.5, color='black')
flierprops = dict(marker='o', markersize=1,
                linestyle='none')
whiskerprops = dict(color='black')
capprops = dict(color='black')
medianprops = dict(linewidth=1.5, linestyle='-', color='black')

bplot1 = ax1.boxplot(bx_plotvals1, labels=names1, notch=False, boxprops=boxprops, whiskerprops=whiskerprops,capprops=capprops, flierprops=flierprops, medianprops=medianprops,showmeans=False)

palette = ['red', 'gray', 'blue', 'yellow','orange']

for xA, valA, c in zip(xs1, vals1, palette):
    ax1.scatter(xA, valA, alpha=0.4, color=c)

sns.set_style("whitegrid")
ind = np.arange(7) 
for ax in fig.get_axes():
    ax.grid(visible=None)
    ax.set_xticks(ticks=ind,labels=['','k5','k7','k10','k15','k20',''],fontsize=15)
    ax.set_ylim(0,2)

ax.set_xlabel('FMH Ksizes')
ax.set_ylabel('dN/dS estimations')
fig.suptitle("dN/dS FMH Across Ksizes",fontsize=20)

fig.savefig(f'{wd}/dNdS_fmh_methods_pairwise.png',bbox_inches='tight')

