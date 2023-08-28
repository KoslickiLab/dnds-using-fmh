#!/usr/bin/python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import rcParams
from scipy import stats
import numpy as np
import math
from sklearn.metrics import mean_absolute_error

def scatter(ground_truth_dnds_results_df,sourmash_compare_dnds_results_df,ksize,wd,gene_length,p_mut):

    str_len=gene_length
    p=p_mut

    expected_dnds_results_df = pd.read_csv(ground_truth_dnds_results_df)
    sourmash_compare_dnds_results_df = pd.read_csv(sourmash_compare_dnds_results_df)

    #column_sourmash_compare_dnds_results_df = 'correct_dnds'
    column_sourmash_compare_dnds_results_df = 'dNdS_ratio'
    column_expected_dnds_results_df = 'dNdS_ground_truth'

    #results = pd.concat([sourmash_compare_dnds_results_df[column_sourmash_compare_dnds_results_df],expected_dnds_results_df[column_expected_dnds_results_df]], axis=1).reindex(expected_dnds_results_df.index)
    results = pd.concat([sourmash_compare_dnds_results_df,expected_dnds_results_df], axis=1).reindex(expected_dnds_results_df.index)
    
    #filter 
    results_2 = results.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio
    results_2 = results_2[results_2['containment_nt']!=results_2['containment_protein']] #remove equal nt and protein containments, this gives HUGE number
    results_2 = results_2[results_2['containment_nt']<=results_2['containment_protein']] #remove if nt containment is greater than protein containment, this gives negative dnds
    Caa=results_2['containment_protein']
    Cnt=results_2['containment_nt']
    k=results_2['ksize']
    results_2 = results_2[Caa-Cnt>=0.001] #remove if too close for p mutation 0.001

    #choosing columns of interest
    results_2 = results_2[[column_sourmash_compare_dnds_results_df,column_expected_dnds_results_df]]
    n=len(results_2[column_expected_dnds_results_df])
    
    """RMSE"""
    MSE=np.square(np.subtract(results_2[column_expected_dnds_results_df],results_2[column_sourmash_compare_dnds_results_df])).mean()
    RMSE=math.sqrt(MSE)

    """MAE"""
    MAE=(abs(results_2[column_expected_dnds_results_df] - results_2[column_sourmash_compare_dnds_results_df]).sum())/len(results_2[column_expected_dnds_results_df])

    """pearson"""
    pcorr  = stats.pearsonr(results_2[column_expected_dnds_results_df],results_2[column_sourmash_compare_dnds_results_df])
    
    """numpy R2 source should tell me how good a model is"""
    corr_matrix=np.corrcoef(results_2[column_expected_dnds_results_df],results_2[column_sourmash_compare_dnds_results_df])
    corr=corr_matrix[0,1]
    R2  = corr**2

    """r2 score from sklearn"""
    #r2_score=r2_score(results_2[column_expected_dnds_results_df],results_2[column_sourmash_compare_dnds_results_df])
    #print(RMSE,corr)
 
    plt.clf()
    plt.scatter(results[column_expected_dnds_results_df],results[column_sourmash_compare_dnds_results_df])

    x_pos=2.2
    plt.text(x_pos,28,f'RMSE={round(RMSE,4)}')
    plt.text(x_pos,26,f'MAE={round(MAE,4)}')
    plt.text(x_pos,24,f'pearson r={round(pcorr[0],4)},p={round(pcorr[1])}')
    plt.text(x_pos,22,f'R^2={round(R2,4)}')
    plt.text(x_pos,20,f'n={n}')
    #plt.text(x_pos,22,f'r2_score={round(r2_score,4)}')

    #plt.text(2.1,28,f'RMSE={RMSE}')
    #plt.text(2.1,26,f'pearson r={corr[0]},p={corr[1]}')

    tmp = [results[column_expected_dnds_results_df].min(), results[column_expected_dnds_results_df].max()]
    plt.plot(tmp, tmp, linestyle='--')
    plt.ylim(-5,30)
    
    #plt.title(f'10,000 nt sequence p=0.1 ksize={ksize} scaled=1')
    plt.title(f'dN/dS real protein-coding sequence scaled=1 p={p} ksize={ksize} len(nt)={str_len}')
    #plt.title(f'dN/dS of K03427 protein-coding genes scaled=1 ksize={ksize} len(nt)={str_len}') 

    #plt.ylabel('smash_containment_dNdS')
    plt.ylabel('smash_containment_dNdS\n outliers removed')
    #plt.ylabel('mahmudur ground truth dN/dD')

    #plt.xlabel('jzr ground truth dN/dS\n(counts total nt differences)')
    plt.xlabel('jzr ground truth dN/dS\n(counts total codon differences)') 

    #plt.savefig(f'{wd}jzr_ground_truth_vs_mahmudur_ground_truth_{p}_{ksize}_{str_len}.png',bbox_inches='tight')
    #plt.savefig(f'{wd}smash_dnds_and_ground_truth_dnds_{p}_{ksize}_{str_len}.png',bbox_inches='tight')
    plt.savefig(f'{wd}smash_dnds_and_ground_truth_dnds_{p}_{ksize}_{str_len}_r2_MAR_RMSE_outliers_rm.png',bbox_inches='tight')
    #plt.savefig(f'{wd}smash_dnds_and_ground_truth_dnds_{p}_{ksize}_{str_len}_r2_MAR_RMSE.png',bbox_inches='tight')

for k in [5,6,7,8,9,10,11,12,13,14,15,20,25,30]:
#for k in [5,6,7,8,9,10,15]:
#for k in [25,30]:

    p_rate=0.1
    length=10002
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=8001
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=4000
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=2000
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=1000
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    p_rate=0.001
    length=10002
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=8001
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=4000
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=2000
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=1000
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    #fmh_dnds=f'{WD}compare_dnds_{k}.csv'
    #scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=fmh_dnds,ksize=k,wd=WD)

    p_rate=0.1
    length=10002
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}_cdn/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=8001
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}_cdn/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=4000
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}_cdn/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=2000
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}_cdn/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=1000
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}_cdn/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    p_rate=0.001
    length=10002
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}_cdn/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=8001
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}_cdn/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=4000
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}_cdn/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=2000
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}_cdn/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    length=1000
    WD=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p_rate}_ksizes_5_30_{length}_cdn/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dnds_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD,gene_length=length,p_mut=p_rate)

    #fmh_dnds=f'{WD}compare_dnds_{k}.csv'
    #scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=fmh_dnds,ksize=k,wd=WD)

