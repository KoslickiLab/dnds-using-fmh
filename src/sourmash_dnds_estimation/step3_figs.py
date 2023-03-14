#!/usr/bin/python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import rcParams

def scatter(ground_truth_dnds_results_df,sourmash_compare_dnds_results_df,ksize,wd):

    expected_dnds_results_df = pd.read_csv(ground_truth_dnds_results_df)
    sourmash_compare_dnds_results_df = pd.read_csv(sourmash_compare_dnds_results_df)

    #results = pd.concat([sourmash_compare_dnds_results_df['dNdS_ratio'],expected_dnds_results_df['dNdS_ground_truth']], axis=1).reindex(expected_dnds_results_df.index)
    
    plt.clf()
    #plt.scatter(results['dNdS_ground_truth'],results['dNdS_ratio'])
    plt.scatter(expected_dnds_results_df['approx_dnds_using_fmh'],sourmash_compare_dnds_results_df['dNdS_ratio'])
    tmp = [min(expected_dnds_results_df['approx_dnds_using_fmh']), max(expected_dnds_results_df['approx_dnds_using_fmh'])]
    plt.plot(tmp, tmp, linestyle='--')
    #plt.ylim(0,20)
    plt.ylim(0,5)
    plt.xlim(0,5)
    plt.title(f'10,000 nt sequence p=0.1 ksize={ksize} scaled=1')
    #plt.title(f'real protein-coding sequence p=0.1 ksize={ksize}')
    #plt.ylabel('smash_containment_dNdS')
    plt.ylabel('jzr_smash_cfrac_dnds')
    plt.xlabel('mahmudur_smash_fmh_dnds')
    #plt.savefig(f'{wd}scatter{ksize}.png',bbox_inches='tight')
    plt.savefig(f'{wd}jzr_cfrac_vs_mahmudur_smash_fmh_dnds{ksize}.png',bbox_inches='tight')

#for k in [2,3,4,5,6,7,8,9,10,11,12,13,14,15]:
#for k in [5,6,7,8,9,10]:
for k in [7]:
    WD='/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/'
    ground_truth=f'{WD}dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'{WD}dNdS_7.csv'
    sourmash_compare_dnds_modified_file=f'{WD}dNdS_7_modified_removed_line_0.csv'
    fmh_dnds=f'{WD}fmh_dnds.txt'
    smash_fmh_dnds=f'{WD}smash_fmh_dnds.txt'
    scatter(ground_truth_dnds_results_df=smash_fmh_dnds,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD)

"""
def multiplot_scatter(ground_truth_dnds_results_df,sourmash_compare_dnds_results_df,ksize):
    expected_dnds_results_df = pd.read_csv(ground_truth_dnds_results_df)
    sourmash_compare_dnds_results_df = pd.read_csv(sourmash_compare_dnds_results_df)

    results = pd.concat([sourmash_compare_dnds_results_df['dNdS_ratio'],expected_dnds_results_df['dNdS_ground_truth']], axis=1).reindex(expected_dnds_results_df.index)
    #minimum = results.min().min()
    #maximum = results.max().max()

    scatter_png = sns.scatterplot(x='dNdS_ratio',y='dNdS_ground_truth',data=results)
    scatter_png = sns.lmplot(data=results, x="dNdS_ratio", y="dNdS_ground_truth")

    plt.title(f'dNdS_ground_truth (koslick_dNdS) vs. Cfrac dNdS ksize={ksize}')
    plt.xlabel('containment_dNdS')

    return(scatter_png)

fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2)
fig.suptitle('A tale of 2 subplots')

ax1.multiplot_scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=5)
ax2.multiplot_scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=6)
ax3.multiplot_scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=7)
ax4.multiplot_scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=8)
ax5.multiplot_scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=9)
ax6.multiplot_scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=10)


WD='/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/'
plt.figure.savefig(f'{WD}test.png')

"""

