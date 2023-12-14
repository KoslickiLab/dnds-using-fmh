import seaborn as sns
import pandas as pd
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

#Distance matrix
labels_df = pd.read_csv(f'{folder1}/compare_protein/compare.prot.{ksize}.mat.labels.txt',sep=',',header=None)
labels_len_list = int(len(list(labels_df[0]))/2)
labels_list = list(labels_df[0])[:labels_len_list]
containments_df = pd.read_csv(f'{folder1}/compare_protein/compare.prot.{ksize}.csv',sep=',')
containments_matrix_df = containments_df[labels_list] .loc[:labels_len_list-1]

containments_matrix_df['species'] = list(containments_matrix_df.columns)
containments_matrix_df = containments_matrix_df.set_index('species',drop=True).rename_axis(index=None)


containments_matrix_df.to_csv('containments_matrix.csv')
print(labels_list)

print(containments_matrix_df.head())


containments_matrix_numpy = containments_matrix_df.to_numpy()
linkage_matrix = linkage(containments_matrix_numpy, "single")

fig, ax = plt.subplots(figsize=(8, 8))

dendrogram=sch.dendrogram(linkage_matrix, labels=labels_list, orientation='right')

fig = plt.gcf()

fig.savefig('dendrogram_test.png',bbox_inches='tight')

