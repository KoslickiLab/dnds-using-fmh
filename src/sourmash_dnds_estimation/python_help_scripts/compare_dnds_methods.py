import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import warnings
from matplotlib.pyplot import figure

#read in pivot tables
methods = ['GNG','GY-HKY','LPB','LWL','NG','YN']
#kaks_calc = pd.read_csv('/data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise_take_2_with_CDS_transcripts/kaks.axt.kaks',sep='\t').pivot_table(index='Sequence',columns='Method')[['Ka/Ks']]
kaks_calc = pd.read_csv('/data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise_take_2_with_CDS_transcripts/kaks.axt.kaks',sep='\t')
kaks_calc['Sequence'] = kaks_calc['Sequence'].replace('>ENSECAT00000016948_Equus.', 'ENSECAT00000016948_Equus. ensembl').replace('>ENSMUST00000111882_Mus.', 'ENSMUST00000111882_Mus. ensembl').replace('>ENSORLT00000022736_Oryzias.', 'ENSORLT00000022736_Oryzias. ensembl').replace('>ENSPPYT00000015088_Pongo.','ENSPPYT00000015088_Pongo. ensembl').replace('>XM_001923765_Danio.','XM_001923765_Danio. refseq').replace('>hsa_7273','hsa:7273 K12567 titin [EC:2.7.11.1] | (RefSeq) TTN, CMD1G, CMH9, CMPD4, CMYP5, EOMFC, HMERF, LGMD2J, LGMDR10, MYLK5, SALMY, TMD; titin (N)')
kaks_calc = kaks_calc.pivot_table(index='Sequence',columns='Method')[['Ka/Ks']]
kaks_calc = kaks_calc['Ka/Ks'][methods].reset_index()

#file is contanation of all dnds_contant files but make sure headers are only in the first line of the file
ksizes = [5,7,10,15,20]
fmh_dnds = pd.read_csv('/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/HIT000324409_pairwise_take_2_with_CDS_transcripts/dnds_constant_all.csv',sep=',').rename(columns={'ksize': 'Method', 'B': 'Sequence'}).pivot_table(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
fmh_dnds = fmh_dnds['dNdS_ratio_constant'][ksizes].reset_index()

#add evola dn/ds modified NG data
evola_dnds = pd.read_csv('/data/jzr5814/sourmash_dnds_estimation/tests/data/evola_data/NFAS/HIT000324409/dnds_HIT000324409.csv',sep='\t')[['Seq.2','dN/dS']].reset_index().rename(columns={'Seq.2': 'Sequence','dN/dS': 'modified_dNdS'})
evola_dnds_method = ['modified_dNdS']

#merge dataframes
df1 = pd.merge(kaks_calc, fmh_dnds, 'left', on = ["Sequence"]).merge(evola_dnds, 'left', on = ["Sequence"])[methods+ksizes+evola_dnds_method]

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

palette = ['red', 'gray', 'blue', 'yellow','orange','purple','brown','pink','olive','cyan','black','darkgreen']

for xA, valA, c in zip(xs1, vals1, palette):
    ax1.scatter(xA, valA, alpha=0.4, color=c)

sns.set_style("whitegrid")
ind = np.arange(14) 
for ax in fig.get_axes():
    ax.grid(visible=None)
    ax.set_xticks(ticks=ind,labels=['','GNG','GY-HKY','LPB','LWL','NG','YN','k5','k7','k10','k15','k20','mod_NG',''],fontsize=15)
    ax.set_ylim(0,1)

ax.set_xlabel('dN/dS methods')
ax.set_ylabel('dN/dS estimations')
fig.suptitle("dN/dS Methods with alignments for KaKs estimations",fontsize=20)
#fig.text(0.07,0.5,'Containment Index',ha='center',va='center',rotation='vertical',fontsize=15)

fig.savefig('/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/HIT000324409_pairwise_take_2_with_CDS_transcripts/dNdS_methods2.png',bbox_inches='tight')

