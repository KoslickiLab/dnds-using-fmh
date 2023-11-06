#!/usr/bin/python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import rcParams
from scipy import stats
import numpy as np
import math
from sklearn.metrics import mean_absolute_error
from scipy.stats import chi2_contingency

def scatter(df1_genomic_cds,df2_genome,ksize,wd):    
    #column_kaks_calc_data = kaks_calc_method
    results = pd.concat([df1_genomic_cds,df2_genome], axis=1)
    
    #filter 
    results = results.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio

    """RMSE"""
    MSE=np.square(np.subtract(results['dnds_genomic_cds'],results['dnds_genomic'])).mean()
    #RMSE=math.sqrt(MSE)

    """MAE"""
    #MAE=(abs(results_2[column_kaks_calc_data] - results_2[column_fracminhash_data]).sum())/len(results_2[column_kaks_calc_data])

    """pearson"""
    pcorr  = stats.pearsonr(results['dnds_genomic_cds'],results['dnds_genomic'])
    
    """numpy R2 source should tell me how good a model is"""
    #corr_matrix=np.corrcoef(results_2[column_kaks_calc_data],results_2[column_fracminhash_data])
    #corr=corr_matrix[0,1]
    #R2  = corr**2

    """chi squared on selectiom"""
    contingency_table = pd.crosstab(results['selection_pressure_genomic_cds'].tolist(),results['selection_pressure_genomic'].tolist())
    chi2, chi2p, dof, expected = chi2_contingency(contingency_table)

    plt.clf()
    plt.scatter(results['dnds_genomic_cds'],results['dnds_genomic'])
    #print(results_2[column_fracminhash_data])
    x_pos=0.1
    plt.text(x_pos,0.22,f'MSE={round(MSE,4)}')
    plt.text(x_pos,0.21,f'pearson r={round(pcorr[0],4)},p={round(pcorr[1])}')
    plt.text(x_pos,0.20,f'chi2={round(chi2,4)}')
    plt.text(x_pos,0.19,f'chi2 p={round(chi2p,4)}')
    #plt.text(x_pos,0.19,f'R^2={round(R2,4)}')
    #plt.text(x_pos,0.17,f'n={n}')
    #plt.text(x_pos,0.23,f'MAE={round(MAE,4)}')
    
    tmp = [results['dnds_genomic_cds'].min(), results['dnds_genomic_cds'].max()]
    plt.plot(tmp, tmp, linestyle='--')
    #plt.ylim(0,2)
    #plt.xlim(0,2)
    
    plt.title(f'dN/dS between Genomic CDS and entire genome scaled=1 ksize={ksize}')
    plt.ylabel('genomic cds dN/dS')
    plt.xlabel(f'entire genome dN/dS')

    plt.savefig(f'{wd}genomic_cds_vs_entire_genome_{ksize}_negative_outliers_rm_2_chi2.png',bbox_inches='tight')

""" Create plot from following data """
#add dnds genomic cds estimations
file1='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_cds_genome/all_dnds_max_containments_with_selection_pressure_fixed_columns_names.csv'
df1 = pd.read_csv(file1,sep=',')[['sequence_comparison','ksize','dNdS_ratio_constant','selection_pressure']].query('0 <= dNdS_ratio_constant <= 2')
df1['sequence_comparison'] = df1['sequence_comparison'].str.replace('_vs_', ',').str.replace('.name_change.fna', '')

#add dnds genomic estimations
file2='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_nt_genome_redo/all_dnds_max_containments_with_selection_pressure_fixed_columns_names.csv'
df2 = pd.read_csv(file2,sep=',')[['sequence_comparison','ksize','dNdS_ratio_constant','selection_pressure']].query('0 <= dNdS_ratio_constant <= 2')
df2['sequence_comparison'] = df2['sequence_comparison'].str.replace('_vs_', ',').str.replace('.name_change.fna', '')

#working directory
WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_nt_genome_redo/'

ksizes = [5]
#ksizes = [7,10]
#ksizes = [15]
#ksizes = [20]
for k in ksizes:
    subset_df1 = df1[df1['ksize'] == k].set_index('sequence_comparison')[['dNdS_ratio_constant','ksize','selection_pressure']].rename(columns={'dNdS_ratio_constant': 'dnds_genomic_cds', 'ksize': 'ksize_genomic_cds','selection_pressure':'selection_pressure_genomic_cds'})
    subset_df2 = df2[df2['ksize'] == k].set_index('sequence_comparison')[['dNdS_ratio_constant','ksize','selection_pressure']].rename(columns={'dNdS_ratio_constant': 'dnds_genomic', 'ksize': 'ksize_genomic','selection_pressure':'selection_pressure_genomic'})
    #FRACMINHASH_DATA=f'{WD}dnds_constant_{k}.csv'
    scatter(df1_genomic_cds=subset_df1,df2_genome=subset_df2,
            ksize=k,wd=WD)
