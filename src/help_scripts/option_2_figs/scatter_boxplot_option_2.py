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
#print(df1)

df1 = df1[(df1 != 0).all(1)]
df1 = df1[df1['0.001-NG'].between(0, 2)]
df1 = df1[df1['0.001-7'].between(0, 2)]
df1 = df1[df1['0.001-10'].between(0, 2)]
df1 = df1[df1['0.001-15'].between(0, 2)]
df1 = df1[df1['0.1-NG'].between(0, 2)]
df1 = df1[df1['0.1-7'].between(0, 2)]
df1 = df1[df1['0.1-10'].between(0, 2)]
df1 = df1[df1['0.1-15'].between(0, 2)]
#print(df1)


####################### Scatter Boxplot #######################

#boxplot
bx_plotvals1, vals1, names1, xs1 = [],[],[],[]
for i, col in enumerate(df1.columns):
    #print(col)
    vals1.append(df1[col].values)
    #print(df1[col].values)
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

palette = ['red', 'gray', 'blue', 'yellow','orange','purple','brown','pink','olive']

for xA, valA, c in zip(xs1, vals1, palette):
    ax1.scatter(xA, valA, alpha=0.4, color=c)

sns.set_style("whitegrid")
ind = np.arange(10) 
for ax in fig.get_axes():
    ax.grid(visible=None)
    ax.set_xticks(ticks=ind,labels=['','0.001-NG','0.1-NG','0.001-7','0.1-7','0.001-10','0.1-10','0.001-15','0.1-15',''],fontsize=15)
    #ax.set_ylim(0,2)

ax.set_xlabel('dN/dS p_rate-method')
ax.set_ylabel('dN/dS estimations')
fig.suptitle("Compare dN/dS Estimates for Option 2",fontsize=20,y=0.91)

#fig.savefig(f'{folder3}/option2_scatter_boxplot.png',bbox_inches='tight')


####################### median values #######################

# Calculate medians for each column
df1 = df1[(df1 > 0) & (df1 < 2)]
medians = df1.median()
means = df1.mean()
maxs = df1.max()
mins = df1.min()

# Create a new DataFrame with the medians
medians_df = pd.DataFrame([medians], index=['Median'])
means_df = pd.DataFrame([means], index=['Mean'])
maxs_df = pd.DataFrame([maxs], index=['Max'])
mins_df = pd.DataFrame([mins], index=['Min'])

medians_df = medians_df.append(means_df)
medians_df = medians_df.append(maxs_df)
medians_df = medians_df.append(mins_df)

# filtering for furtehr analysis
above1 = df1[(df1 > 1) & (df1 < 2)].replace([np.inf, -np.inf], np.nan).dropna()
below1 = df1[(df1 > 0) & (df1 < 1)].replace([np.inf, -np.inf], np.nan).dropna()

# Calculate medians for each column
above1_medians = above1.median()
above1_means = above1.mean()
above1_maxs = above1.max()
above1_mins = above1.min()

above1_medians_df = pd.DataFrame([above1_medians], index=['above1_medians'])
above1_means_df = pd.DataFrame([above1_means], index=['above1_means'])
above1_maxs_df = pd.DataFrame([above1_maxs], index=['above1_maxs'])
above1_mins_df = pd.DataFrame([above1_mins], index=['above1_mins'])

below1_medians = below1.median()
below1_means = below1.mean()
below1_maxs = below1.max()
below1_mins = below1.min()

below1_medians_df = pd.DataFrame([below1_medians], index=['below1_medians'])
below1_means_df = pd.DataFrame([below1_means], index=['below1_means'])
below1_maxs_df = pd.DataFrame([below1_maxs], index=['below1_maxs'])
below1_mins_df = pd.DataFrame([below1_mins], index=['below1_mins'])

medians_df = medians_df.append(above1_medians_df)
medians_df = medians_df.append(above1_means_df)
medians_df = medians_df.append(above1_maxs_df)
medians_df = medians_df.append(above1_mins_df)

medians_df = medians_df.append(below1_medians_df)
medians_df = medians_df.append(below1_means_df)
medians_df = medians_df.append(below1_maxs_df)
medians_df = medians_df.append(below1_mins_df)

medians_df.round(5).to_csv(f'{folder3}/stats.csv')