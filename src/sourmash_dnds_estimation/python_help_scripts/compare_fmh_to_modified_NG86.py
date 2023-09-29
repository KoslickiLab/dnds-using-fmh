#!/usr/bin/python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import rcParams
from scipy import stats
import numpy as np
import math
from sklearn.metrics import mean_absolute_error

def scatter(ground_truth_data,fracminhash_data,ksize,wd,kaks_calc_method):

    #read modified NG
    kaks_calc_data = ground_truth_data[['Sequence',kaks_calc_method]]
    #read in fracminhash dataframe
    fracminhash_data = pd.read_csv(fracminhash_data)

    #choose column names depending on the data
    column_fracminhash_data = 'dNdS_ratio_constant' #USING SOURMASH
    
    column_kaks_calc_data = kaks_calc_method

    results = pd.concat([fracminhash_data,kaks_calc_data], axis=1).reindex(kaks_calc_data.index)
    
    #filter 
    results_2 = results.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio

    #choosing columns of interest
    results_2 = results_2[[column_fracminhash_data,column_kaks_calc_data]]
    results_2 = results_2[(results_2[column_fracminhash_data] >= 0)]
    print(results_2)
    n=len(results_2[column_kaks_calc_data])
    
    """RMSE"""
    MSE=np.square(np.subtract(results_2[column_kaks_calc_data],results_2[column_fracminhash_data])).mean()
    RMSE=math.sqrt(MSE)

    """MAE"""
    MAE=(abs(results_2[column_kaks_calc_data] - results_2[column_fracminhash_data]).sum())/len(results_2[column_kaks_calc_data])

    """pearson"""
    pcorr  = stats.pearsonr(results_2[column_kaks_calc_data],results_2[column_fracminhash_data])
    
    """numpy R2 source should tell me how good a model is"""
    corr_matrix=np.corrcoef(results_2[column_kaks_calc_data],results_2[column_fracminhash_data])
    corr=corr_matrix[0,1]
    R2  = corr**2
 
    plt.clf()
    plt.scatter(results_2[column_kaks_calc_data],results_2[column_fracminhash_data])
    print(results_2[column_fracminhash_data])
    x_pos=0.3
    plt.text(x_pos,0.25,f'RMSE={round(RMSE,4)}')
    plt.text(x_pos,0.23,f'MAE={round(MAE,4)}')
    plt.text(x_pos,0.21,f'pearson r={round(pcorr[0],4)},p={round(pcorr[1])}')
    plt.text(x_pos,0.19,f'R^2={round(R2,4)}')
    plt.text(x_pos,0.17,f'n={n}')

    tmp = [results_2[column_kaks_calc_data].min(), results_2[column_kaks_calc_data].max()]
    plt.plot(tmp, tmp, linestyle='--')
    #plt.ylim(0,5)
    
    plt.title(f'dN/dS between real protein genes scaled=1 ksize={ksize}')
    plt.ylabel('fmh dN/dS')
    plt.xlabel(f'{kaks_calc_method} dN/dS (EVOLA data)')

    plt.savefig(f'{wd}fmh_constant_and_{kaks_calc_method}_{ksize}_negative_outliers_rm_2.png',bbox_inches='tight')

""" Create plot from following data """
#add fmh dn/ds data
ksizes = [5,7,10,15,20]
fmh_dnds = pd.read_csv('/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/HIT000324409/dnds_constant_all.csv',sep=',').rename(columns={'ksize': 'Method', 'B': 'Sequence'}).pivot(index='Sequence',columns='Method')[['dNdS_ratio_constant']]
fmh_dnds = fmh_dnds['dNdS_ratio_constant'][ksizes].reset_index()
#add evola dn/ds modified NG data
evola_dnds = pd.read_csv('/data/jzr5814/sourmash_dnds_estimation/tests/data/evola_data/NFAS/HIT000324409/dnds_HIT000324409.csv',sep='\t')[['Seq.2','dN/dS']].reset_index().rename(columns={'Seq.2': 'Sequence','dN/dS': 'modified_NG'})
evola_dnds_method = ['modified_NG']
#merge dataframes
df1 = pd.merge(evola_dnds, fmh_dnds, 'left', on = ["Sequence"])[['Sequence']+ksizes+evola_dnds_method]
#working directory
WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/HIT000324409/'

for k in ksizes:
    FRACMINHASH_DATA=f'{WD}dnds_constant_{k}.csv'
    scatter(ground_truth_data=df1,fracminhash_data=FRACMINHASH_DATA,
            ksize=k,wd=WD,
            kaks_calc_method='modified_NG')
