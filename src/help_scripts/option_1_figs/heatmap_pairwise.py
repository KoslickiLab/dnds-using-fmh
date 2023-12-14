"""
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Example data (replace this with your actual data)
data = pd.DataFrame({
    'Category': ['A', 'B', 'C'],
    'Column1': [1, 2, 3],
    'Column2': [4, 5, 6],
    'Column3': [7, 8, 9]
})

print(data)

"""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pylab import savefig
import numpy as np

import pickle
import warnings

import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import linkage, dendrogram 

import matplotlib.pyplot as plt



title='Compare Selection Pressures for Option 2'
folder1="/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/marine_bacteria"
fmhdnds='all_dnds_constant.csv'

ksize=15

#dN/dS values
fmh_dnds = pd.read_csv(f'{folder1}/{fmhdnds}',sep=',')[['A','B','ksize','containment_nt','containment_protein','dNdS_ratio_constant']]
sorted_df = fmh_dnds[fmh_dnds['ksize'] == ksize].sort_values(by=['A', 'B'])
print(sorted_df)

dic={}
for index, row in sorted_df.iterrows():
    gene1_value = row['A']
    col3_value = row['dNdS_ratio_constant']
    if gene1_value not in dic:
            dic[gene1_value] = []
            dic[gene1_value].append(col3_value)
    elif gene1_value in dic:
            dic[gene1_value].append(col3_value)
dNdS_df = pd.DataFrame(dic).fillna(1)
dNdS_df['species'] = list(dNdS_df.columns)
dNdS_df = dNdS_df.set_index('species',drop=True).rename_axis(index=None)

plt=sns.heatmap(dNdS_df,cmap='seismic')
plt.figure.savefig('dnds_test.png',bbox_inches='tight')

dNdS_df.to_csv('dNdS_matrix.csv')
