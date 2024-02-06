#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse, time, os, subprocess
from dNdS import reportCI,CfracdNdS
import pandas as pd
pd.set_option('display.max_columns', None)

def main(args):

    ###Arguments
    WD = args.wd
    nt_cfracs = args.nt
    protein_cfracs = args.protein
    ksize = args.k

    # for both nt and protein containments, make sure to run: awk 'NR == 1 || FNR > 1' nt_containment* > nt_all_cfrac.csv
    nt_df = pd.read_csv(nt_cfracs,sep=',')
    print(nt_df.describe())
    prot_df = pd.read_csv(protein_cfracs,sep=',')

    dna_cfrac_median_list = []
    ANI_approx_median_list = []
    prot_cfrac_median_list = []
    AAI_approx_median_list = []

    ### obtain median cfrac per ksize
    klist = ksize.split(',')
    print(klist)
    for k in klist:
        k=int(k)
        nt_df_temp = nt_df[nt_df['ksize'] == k]
        dna_cfrac_median_list.append(nt_df_temp['DNA_Cfrac'].median())
        ANI_approx_median_list.append(nt_df_temp['approx_ANI'].median())

        prot_df_temp = prot_df[prot_df['ksize'] == k]
        prot_cfrac_median_list.append(prot_df_temp['protein_Cfrac'].median())
        AAI_approx_median_list.append(prot_df_temp['approx_AAI'].median())

    #produce dictionary for dataframe creation
    median_dict={'ksize':klist,'approx_ANI_median':ANI_approx_median_list,'DNA_Cfrac_median':dna_cfrac_median_list,'approx_AAI_median':AAI_approx_median_list,'protein_Cfrac_median':prot_cfrac_median_list}

    #report as csv file    
    median_df = pd.DataFrame.from_dict(median_dict)
    median_df.to_csv(f'{WD}cfracs_and_approximations_median.csv')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Report medians for Cfracs and AXI approximations'
    )

    parser.add_argument(
        '--k',
        type=str,
        help = 'Identify a list of ksizes.\
        ksize is required to be the same used in both the containment indexes calculated for nucleotide and protein sequences.'
    )

    parser.add_argument(
        '--nt',
        help = 'Input of nucleotide compare.csv file that has been modifed and includes all ksizes for analysis.\
        This file is a pairwise matrix produced from sourmash compare that includes containment indexes between nucleotide sequences.'
    )

    parser.add_argument(
        '--protein',
        help = 'Input of protein compare.csv file that has been modifed and includes all ksizes for analysis.\
        This file is a pairwise matrix produced from sourmash compare that includes containment indexes between protein sequences.'
    )

    parser.add_argument(
    '--wd',
    help = 'Output directory for CSV file with dN/dS ratio estimates between sequences evaluated.'
    )

    args = parser.parse_args()

    main(args)
