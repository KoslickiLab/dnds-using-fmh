import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pylab import savefig
import numpy as np

title='Max C(A,B) for 4,000 nt long sequence with p rate 0.001 Compared to Other dN/dS Methods'
folder1="real_data_0.001_4000_cdn"
folder2="real_data_0.001_ksizes_5_30_4000_cdn"


#read in pivot tables
#methods = ['GNG','GY-HKY','LPB','LWL','NG','YN']
methods = ['GY-HKY','YN','GNG','NG']


kaks_calc = pd.read_csv(f'/data/jzr5814/kaks_calc_tool_analysis/{folder1}/kaks_sequences.axt.kaks',sep='\t')
#kaks_calc['Sequence'] = kaks_calc['Sequence'].replace('>ENSECAT00000016948_Equus.', 'ENSECAT00000016948_Equus. ensembl').replace('>ENSMUST00000111882_Mus.', 'ENSMUST00000111882_Mus. ensembl').replace('>ENSORLT00000022736_Oryzias.', 'ENSORLT00000022736_Oryzias. ensembl').replace('>ENSPPYT00000015088_Pongo.','ENSPPYT00000015088_Pongo. ensembl').replace('>XM_001923765_Danio.','XM_001923765_Danio. refseq').replace('>hsa_7273','hsa:7273 K12567 titin [EC:2.7.11.1] | (RefSeq) TTN, CMD1G, CMH9, CMPD4, CMYP5, EOMFC, HMERF, LGMD2J, LGMDR10, MYLK5, SALMY, TMD; titin (N)')
#kaks_calc = kaks_calc.pivot_table(index='Sequence',columns='Method')[['Ka/Ks']]
kaks_calc = kaks_calc.pivot(index='Sequence',columns='Method')[['Ka/Ks']] #use when no duplicates
#kaks_calc = kaks_calc['Ka/Ks'][methods].reset_index()
kaks_calc = kaks_calc['Ka/Ks']

#file is contanation of all dnds_contant files but make sure headers are only in the first line of the file
ksizes = [5,7,10,15,20]
fmh_dnds = pd.read_csv(f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/{folder2}/all_dnds_max_containments.csv',sep=',').rename(columns={'ksize': 'Method', 'B': 'Sequence'}).pivot_table(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
#fmh_dnds = pd.read_csv(f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/{folder}/all_dnds_max_containments.csv',sep=',').rename(columns={'ksize': 'Method', 'sequence_comparison': 'Sequence'}).pivot(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
#fmh_dnds = fmh_dnds['dNdS_ratio_constant'][ksizes].reset_index()
fmh_dnds = fmh_dnds['dNdS_ratio_constant'][ksizes]


#merge dataframes
df1 = pd.merge(kaks_calc, fmh_dnds, 'left', on = ["Sequence"])
df1.to_csv(f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/{folder2}/dNdS_methods_pairwise_without_mean_without_evola_maxC.csv')
print(df1)

df1=df1[methods+ksizes]

#filter outliers
df1 = df1.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio
df1 = df1[(df1 != 0).all(1)]

plt = sns.heatmap(df1, vmin=0, vmax=2,xticklabels=True, cmap='seismic')
plt.set_title(title)
plt.figure.savefig(f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/{folder2}/heatmap_pairwise_dNdS_maxC.png',bbox_inches='tight')
