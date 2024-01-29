#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse, time, os, subprocess
from dNdS import reportCI,CfracdNdS
import pandas as pd
pd.set_option('display.max_columns', None)

def main(args):

    ###Arguments
    WD = args.wd
    nt_compare = args.nt
    protein_compare = args.protein
    ksize = args.k

    ### Obtain containment from matrix files produced by sourmash compare
    nt_df = reportCI.grab_containment_from_mat_ground_truth(mat_df=nt_compare,ksize=ksize)
    nt_df = nt_df[nt_df['A'] == 'ref_gene'].rename(columns={"containment": "DNA_Cfrac"})
    nt_df['approx_ANI'] = nt_df['DNA_Cfrac']**(1/nt_df['ksize'])
    print(nt_df)

    protein_df = reportCI.grab_containment_from_mat_ground_truth(mat_df=protein_compare,ksize=ksize)
    protein_df = protein_df[protein_df['A'] == 'ref_gene'].rename(columns={"containment": "protein_Cfrac"})
    protein_df['approx_AAI'] = protein_df['protein_Cfrac']**(1/protein_df['ksize'])

    nt_df.to_csv(f'{WD}ANI_approximations_{ksize}.csv')
    protein_df.to_csv(f'{WD}AAI_approximations_{ksize}.csv')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'dN/dS estimator for metagenomic data using the containment index between k-mer sets of genomic samples'
    )

    parser.add_argument(
        '--k',
        type=int,
        help = 'Identify the ksize used to produce containment indexes between two sequences.\
        ksize is required to be the same used in both the containment indexes calculated for nucleotide and protein sequences.'
    )

    parser.add_argument(
        '--nt',
        help = 'Input of nucleotide compare.csv file.\
        This file is a pairwise matrix produced from sourmash compare that includes containment indexes between nucleotide sequences.'
    )

    parser.add_argument(
        '--protein',
        help = 'Input of protein compare.csv file.\
        This file is a pairwise matrix produced from sourmash compare that includes containment indexes between protein sequences.'
    )

    parser.add_argument(
    '--o',
    help = 'Output CSV file with dN/dS ratio estimates between sequences evaluated. The CSV file reports sequence A, sequence B,\
    ksize, and containment index'
    )

    parser.add_argument(
    '--wd',
    help = 'Output directory for CSV file with dN/dS ratio estimates between sequences evaluated.'
    )

    args = parser.parse_args()

    main(args)
