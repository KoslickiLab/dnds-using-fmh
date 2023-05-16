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
    kaks_calc_data = kaks_calc_data[kaks_calc_data['Method']==kaks_calc_method].tail(-1).reset_index()[['Sequence','Ka/Ks']]

    #read in fracminhash dataframe
    fracminhash_data = pd.read_csv(fracminhash_data)

    #choose column names depending on the data
    #column_fracminhash_data = 'dNdS_ground_truth' #USING PYTHON GROUND TRUTH
    column_fracminhash_data = 'dNdS_ratio' #USING SOURMASH
    
    column_kaks_calc_data = 'Ka/Ks'

    #results = pd.concat([fracminhash_data,kaks_calc_data], axis=1).reindex(kaks_calc_data.index)
    results = pd.concat([fracminhash_data,kaks_calc_data], axis=1).reindex(kaks_calc_data.index)
    #results.to_csv(f"{WD}raw_compare_groundtruthfracminhash_KaKs.csv")
    #results.to_csv(f"{WD}raw_compare_groundtruthfracminhash_KaKs{kaks_calc_method}.csv")
    
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
    MSE=np.square(np.subtract(results_2[column_kaks_calc_data],results_2[column_fracminhash_data]*(0.77/2.23))).mean()
    RMSE=math.sqrt(MSE)

    """MAE"""
    MAE=(abs(results_2[column_kaks_calc_data] - results_2[column_fracminhash_data]*(0.77/2.23)).sum())/len(results_2[column_kaks_calc_data])

    """pearson"""
    pcorr  = stats.pearsonr(results_2[column_kaks_calc_data],results_2[column_fracminhash_data]*(0.77/2.23))
    
    """numpy R2 source should tell me how good a model is"""
    corr_matrix=np.corrcoef(results_2[column_kaks_calc_data],results_2[column_fracminhash_data]*(0.77/2.23))
    corr=corr_matrix[0,1]
    R2  = corr**2

    """r2 score from sklearn"""
    #r2_score=r2_score(results_2[column_kaks_calc_data],results_2[column_fracminhash_data])
    #print(RMSE,corr)
 
    plt.clf()
    plt.scatter(results[column_kaks_calc_data],results[column_fracminhash_data]*(0.77/2.23))

    x_pos=0.25
    plt.text(x_pos,56,f'RMSE={round(RMSE,4)}')
    plt.text(x_pos,53,f'MAE={round(MAE,4)}')
    plt.text(x_pos,50,f'pearson r={round(pcorr[0],4)},p={round(pcorr[1])}')
    plt.text(x_pos,47,f'R^2={round(R2,4)}')
    plt.text(x_pos,44,f'n={n}')
    #plt.text(x_pos,22,f'r2_score={round(r2_score,4)}')

    #plt.text(2.1,28,f'RMSE={RMSE}')
    #plt.text(2.1,26,f'pearson r={corr[0]},p={corr[1]}')

    tmp = [results[column_kaks_calc_data].min(), results[column_kaks_calc_data].max()]
    plt.plot(tmp, tmp, linestyle='--')
    plt.ylim(-5,60)
    
    #plt.title(f'10,000 nt sequence p={p_mut} ksize={ksize} scaled=1')
    plt.title(f'dN/dS real protein-coding sequence scaled=1 p={p_mut} ksize={ksize} len(nt)={gene_length}')
    #plt.title(f'dN/dS of K03427 protein-coding genes scaled=1 ksize={ksize} len(nt)={str_len}') 

    plt.ylabel('fmh dN/dS\nsourmash containment')
    #plt.ylabel('jzr ground truth dN/dS\n(counts total nt differences)')
    #plt.ylabel('jzr ground truth dN/dS\n(counts total codon differences)')

    plt.xlabel(f'{kaks_calc_method} dN/dS\nKa/Ks calculator')

    #plt.savefig(f'{wd}{kaks_calc_method}_vs_jzr_ground_truth_{p_mut}_{ksize}_{gene_length}.png',bbox_inches='tight')
    plt.savefig(f'{wd}fmh_constant_and_{kaks_calc_method}_{p_mut}_{ksize}_{gene_length}_r2_MAE_RMSE_outliers_rm.png',bbox_inches='tight')

for k in [5,6,7,8,9,10,11,12,13,14,15,20,25,30]:

    p_rate=0.001
    length=4000
    #KAKS_CALC_METHOD='GY-HKY'
    KAKS_CALC_METHOD='NG'
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}_cdn/'
    GROUND_TRUTH_FRACMINHASH=f'{WD}dNdS_ground_truth.csv' #compare to jzr ground truth
    FRACMINHASH_DATA=f'{WD}dnds_{k}.csv'
    ka_ks_data = f'/data/jzr5814/kaks_calc_tool_analysis/real_data_{p_rate}_{length}_cdn/kaks_sequences.axt.kaks'

    scatter(ground_truth_data=ka_ks_data,fracminhash_data=FRACMINHASH_DATA,
            ksize=k,wd=WD,gene_length=length,p_mut=p_rate,
            kaks_calc_method=KAKS_CALC_METHOD)
