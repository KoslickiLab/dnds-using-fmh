#!/usr/bin/python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import rcParams
from scipy import stats
import numpy as np
import math
from sklearn.metrics import mean_absolute_error

def scatter(ground_truth_data,fracminhash_data,ksize,wd,gene_length,p_mut,kaks_calc_method):

    #fix kaks calculator dataframe by choosing method of interest, removing refvsref dnds and reindexing
    kaks_calc_data = pd.read_csv(ground_truth_data,delimiter='\t')
    kaks_calc_data = kaks_calc_data[kaks_calc_data['Method']==kaks_calc_method].tail(-1).reset_index()[['Sequence','Ka']]

    #read in fracminhash dataframe
    fracminhash_data = pd.read_csv(fracminhash_data)

    #choose column names depending on the data
    column_fracminhash_data = 'dN_constant' #USING SOURMASH
    
    column_kaks_calc_data = 'Ka'

    results = pd.concat([fracminhash_data,kaks_calc_data], axis=1).reindex(kaks_calc_data.index)
    
    #filter 
    results_2 = results.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio
    
    #filter when using sourmash results
    results_2 = results_2[results_2['containment_nt']!=results_2['containment_protein']] #remove equal nt and protein containments, this gives HUGE number
    results_2 = results_2[results_2['containment_nt']<=results_2['containment_protein']] #remove if nt containment is greater than protein containment, this gives negative dnds
    Caa=results_2['containment_protein']
    Cnt=results_2['containment_nt']
    results_2 = results_2[Caa-Cnt>=0.001] #remove if too close for p mutation 0.001
    
    #choosing columns of interest
    results_2 = results_2[[column_fracminhash_data,column_kaks_calc_data]]
    
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
    plt.scatter(results[column_kaks_calc_data],results[column_fracminhash_data])
    """
    x_pos=0.0004
    plt.text(x_pos,0.0023,f'RMSE={round(RMSE,4)}')
    plt.text(x_pos,0.0022,f'MAE={round(MAE,4)}')
    plt.text(x_pos,0.0021,f'pearson r={round(pcorr[0],4)},p={round(pcorr[1])}')
    plt.text(x_pos,0.0020,f'R^2={round(R2,4)}')
    plt.text(x_pos,0.0019,f'n={n}')
    """
    tmp = [results[column_kaks_calc_data].min(), results[column_kaks_calc_data].max()]
    plt.plot(tmp, tmp, linestyle='--')
    #plt.ylim(0,0.0025)
    
    plt.title(f'dN real protein-coding sequence scaled=1 p={p_mut} ksize={ksize} len(nt)={gene_length}')
    plt.ylabel('fmh dN\nsourmash containment')
    plt.xlabel(f'{kaks_calc_method} dN\nKa/Ks calculator')

    plt.savefig(f'{wd}Ka_fmh_constant_and_{kaks_calc_method}_{p_mut}_{ksize}_{gene_length}_r2_MAE_RMSE_outliers_rm.png',bbox_inches='tight')

for k in [5,6,7,8,9,10,11,12,13,14,15,20,25,30]:

    p_rate=0.1
    length=10002
    #KAKS_CALC_METHOD='GY-HKY'
    KAKS_CALC_METHOD='NG'
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}_cdn/'
    FRACMINHASH_DATA=f'{WD}dnds_constant_{k}.csv'
    ka_ks_data = f'/data/jzr5814/kaks_calc_tool_analysis/real_data_{p_rate}_{length}_cdn/kaks_sequences.axt.kaks'

    scatter(ground_truth_data=ka_ks_data,fracminhash_data=FRACMINHASH_DATA,
            ksize=k,wd=WD,gene_length=length,p_mut=p_rate,
            kaks_calc_method=KAKS_CALC_METHOD)
