#!/usr/bin/python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import rcParams
from scipy import stats
import numpy as np
import math

def scatter(ground_truth_dnds_results_df,sourmash_compare_dnds_results_df,ksize,wd):

    str_len=1000
    p=0.1

    expected_dnds_results_df = pd.read_csv(ground_truth_dnds_results_df)
    sourmash_compare_dnds_results_df = pd.read_csv(sourmash_compare_dnds_results_df)

    column_sourmash_compare_dnds_results_df = 'correct_dnds'
    #column_sourmash_compare_dnds_results_df = 'dNdS_ratio'
    column_expected_dnds_results_df = 'dNdS_ground_truth'

    results = pd.concat([sourmash_compare_dnds_results_df[column_sourmash_compare_dnds_results_df],expected_dnds_results_df[column_expected_dnds_results_df]], axis=1).reindex(expected_dnds_results_df.index)
    results_2 = results.replace([np.inf, -np.inf], np.nan).dropna()
    #print(results_2)
    #print(results_2.query('correct_dnds != dNdS_ground_truth'))

    MSE=np.square(np.subtract(results_2[column_expected_dnds_results_df],results_2[column_sourmash_compare_dnds_results_df])).mean()
    RMSE=math.sqrt(MSE)
    corr  = stats.pearsonr(results_2[column_expected_dnds_results_df],results_2[column_sourmash_compare_dnds_results_df])
    #print(RMSE,corr)
 
    plt.clf()
    plt.scatter(results[column_expected_dnds_results_df],results[column_sourmash_compare_dnds_results_df])

    #plt.text(0,25,f'RMSE={round(RMSE,4)}')
    #plt.text(0,24,f'pearson r={round(corr[0],4)},p={round(corr[1])}')

    plt.text(2,4,f'RMSE={RMSE}')
    plt.text(2,3.9,f'pearson r={corr[0]},p={corr[1]}')

    tmp = [results[column_expected_dnds_results_df].min(), results[column_expected_dnds_results_df].max()]
    plt.plot(tmp, tmp, linestyle='--')
    #plt.ylim(-5,30)
    
    #plt.title(f'10,000 nt sequence p=0.1 ksize={ksize} scaled=1')
    plt.title(f'dN/dS real protein-coding sequence scaled=1 p={p} ksize={ksize} len(nt)={str_len}')
    #plt.title(f'dN/dS of K03427 protein-coding genes scaled=1 ksize={ksize} len(nt)={str_len}') 

    #plt.ylabel('smash_containment_dNdS')
    plt.ylabel('mahmudur ground truth dN/dD')

    plt.xlabel('jzr ground truth dN/dS\n(counts total nt differences)')
    #plt.xlabel('jzr ground truth dN/dS\n(counts total codon differences)') 

    plt.savefig(f'{wd}jzr_ground_truth_vs_mahmudur_ground_truth_{p}_{ksize}_{str_len}.png',bbox_inches='tight')
    #plt.savefig(f'{wd}smash_dnds_and_ground_truth_dnds_{p}_{ksize}_{str_len}.png',bbox_inches='tight')

for k in [5,6,7,8,9,10,11,12,13,14,15,20,25,30]:
#for k in [5,6,7,8,9,10,15]:
#for k in [25,30]:
    WD='/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_0.1_ksizes_5_30_1000/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'

    #sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    #scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD)

    fmh_dnds=f'{WD}compare_dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=fmh_dnds,ksize=k,wd=WD)

