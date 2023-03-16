#!/usr/bin/python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import rcParams

for ksize in [2,3,4,5,6,7,8,9,10,11,12,13,14,15,20]:
    WD='/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/random_10000_nt_sequence_compared_to_mahmudur_code_0.1_ksizes_2_to_20/'
    containment_comparisons=f'{WD}containment_comparisons{ksize}.csv'
    containment_comparisons_df = pd.read_csv(containment_comparisons)

    plt.clf()
    plt.scatter(containment_comparisons_df['mahmudur_p_nt'],containment_comparisons_df['jzr_p_nt'])
    tmp = [min(containment_comparisons_df['mahmudur_p_nt']), max(containment_comparisons_df['mahmudur_p_nt'])]
    plt.plot(tmp, tmp, linestyle='--')
    #plt.ylim(0,20)
    #plt.ylim(0,5)
    #plt.xlim(0,5)
    plt.title(f'nt containment 10,000 nt sequence p=0.1 ksize={ksize} scaled=1')
    #plt.title(f'real protein-coding sequence p=0.1 ksize={ksize}')
    #plt.ylabel('smash_containment_dNdS')
    plt.ylabel('jzr using sourmash')
    plt.xlabel('mahmudur using python')
    plt.savefig(f'{WD}jzr_cfrac_vs_mahmudur_nt_containment{ksize}.png',bbox_inches='tight')