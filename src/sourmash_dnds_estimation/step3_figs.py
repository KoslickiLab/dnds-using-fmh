#!/usr/bin/python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import rcParams
from scipy import stats
import numpy as np
import math

def scatter(ground_truth_dnds_results_df,sourmash_compare_dnds_results_df,ksize,wd):

    str_len=10002
    p=0.1

    expected_dnds_results_df = pd.read_csv(ground_truth_dnds_results_df)
    sourmash_compare_dnds_results_df = pd.read_csv(sourmash_compare_dnds_results_df)

    column_sourmash_compare_dnds_results_df = 'correct_dnds'
    column_expected_dnds_results_df = 'dNdS_ground_truth'

    results = pd.concat([sourmash_compare_dnds_results_df[column_sourmash_compare_dnds_results_df],expected_dnds_results_df[column_expected_dnds_results_df]], axis=1).reindex(expected_dnds_results_df.index)
    results_2 = results.replace([np.inf, -np.inf], np.nan).dropna()
    print(results_2)

    MSE=np.square(np.subtract(results_2[column_expected_dnds_results_df],results_2[column_sourmash_compare_dnds_results_df])).mean()
    RMSE=math.sqrt(MSE)
    corr  = stats.pearsonr(results_2[column_expected_dnds_results_df],results_2[column_sourmash_compare_dnds_results_df])
    print(RMSE,corr)
 
    plt.clf()
    plt.scatter(results[column_expected_dnds_results_df],results[column_sourmash_compare_dnds_results_df])

    plt.text(2.2,2.9,f'RMSE={RMSE}')
    #plt.text(2,5.3,f'pearson r={round(corr[0],4)},p={round(corr[1])}')
    plt.text(2.2,2.87,f'pearson r={corr[0]},p={corr[1]}')

    #plt.scatter(expected_dnds_results_df['approx_dnds_using_fmh'],sourmash_compare_dnds_results_df['dNdS_ratio'])
    #plt.scatter(expected_dnds_results_df['dNdS_ground_truth'],sourmash_compare_dnds_results_df['dNdS_ratio'])
    tmp = [results[column_expected_dnds_results_df].min(), results[column_expected_dnds_results_df].max()]
    plt.plot(tmp, tmp, linestyle='--')
    #plt.ylim(-125,125)
    
    #plt.title(f'10,000 nt sequence p=0.1 ksize={ksize} scaled=1')
    plt.title(f'dN/dS real protein-coding sequence scaled=1 p={p} ksize={ksize} len(nt)={str_len}')
    #plt.title(f'dN/dS of K03427 protein-coding genes scaled=1 ksize={ksize} len(nt)={str_len}') 
    #plt.ylabel('smash_containment_dNdS')
    plt.ylabel('mahmudur groun truth dN/dD')
    plt.xlabel('jzr ground truth dN/dS') 
    plt.savefig(f'{wd}jzr_ground_truth_vs_mahmudur_ground_truth_{p}_{ksize}_{str_len}.png',bbox_inches='tight')
    #plt.savefig(f'{wd}smash_dnds_and_ground_truth_dnds_{p}_{ksize}_{str_len}.png',bbox_inches='tight'))

#for k in [5,6,7,8,25,30]:
for k in [5,6,7,8,9,10,15,20,25,30]:
#for k in [20,25]:
    WD='/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_0.1_ksizes_5_30_10002_cdn_differences/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    #sourmash_compare_dnds=f'{WD}dNdS_{k}.csv'
    #sourmash_compare_dnds_modified_file=f'{WD}dNdS_7_modified_removed_line_0.csv'
    fmh_dnds=f'{WD}compare_dnds_{k}.csv'
    #smash_fmh_dnds=f'{WD}smash_fmh_dnds.txt'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=fmh_dnds,ksize=k,wd=WD)


