#!/usr/bin/python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import rcParams

def scatter(ground_truth_dnds_results_df,sourmash_compare_dnds_results_df,p,ksize,wd):

    expected_dnds_results_df = pd.read_csv(ground_truth_dnds_results_df)
    sourmash_compare_dnds_results_df = pd.read_csv(sourmash_compare_dnds_results_df)

    results = pd.concat([sourmash_compare_dnds_results_df['dNdS_ratio'],expected_dnds_results_df['dNdS_ground_truth']], axis=1).reindex(expected_dnds_results_df.index)
    
    plt.clf()
    plt.scatter(results['dNdS_ground_truth'],results['dNdS_ratio'])
    #plt.scatter(expected_dnds_results_df['approx_dnds_using_fmh'],sourmash_compare_dnds_results_df['dNdS_ratio'])
    #plt.scatter(expected_dnds_results_df['dNdS_ground_truth'],sourmash_compare_dnds_results_df['dNdS_ratio'])
    tmp = [min(results['dNdS_ground_truth']), max(results['dNdS_ground_truth'])]
    #tmp = [min(expected_dnds_results_df['approx_dnds_using_fmh']), max(expected_dnds_results_df['approx_dnds_using_fmh'])]
    plt.plot(tmp, tmp, linestyle='--')
    #plt.ylim(0,20)
    #plt.ylim(0,5)
    #plt.xlim(0,5)
    #plt.title(f'10,000 nt sequence p=0.1 ksize={ksize} scaled=1')
    plt.title(f'dN/dS real protein-coding sequence (1,470 nt) scaled=1 p={p} ksize={ksize}')
    #plt.ylabel('smash_containment_dNdS')
    plt.ylabel('sourmash dN/dD')
    plt.xlabel('ground truth dN/dS')
    #plt.savefig(f'{wd}scatter{ksize}.png',bbox_inches='tight')
    plt.savefig(f'{wd}smash_dnds_and_ground_truth_dnds{ksize}.png',bbox_inches='tight')

for k in [5,6,7,8,9,10]:
#for k in [7]:
    WD='/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_0.1_ksizes_5_10/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dNdS_{k}.csv'
    #sourmash_compare_dnds_modified_file=f'{WD}dNdS_7_modified_removed_line_0.csv'
    #fmh_dnds=f'{WD}compare_dnds_{k}.csv'
    #smash_fmh_dnds=f'{WD}smash_fmh_dnds.txt'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,p=0.1,ksize=k,wd=WD)


