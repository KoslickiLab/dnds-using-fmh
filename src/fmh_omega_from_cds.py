#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse, time, os, subprocess
from fmh_omega import helperfuncs,dnds
import pandas as pd
pd.set_option('display.max_columns', None)

def main(args):

    ###Arguments
    WD = args.wd
    nt_compare = args.nt
    protein_compare = args.protein
    ksize = args.k
    output = args.o

    ### Obtain containment from matrix files produced by sourmash compare
    nt_df = helperfuncs.containments(mat_df=nt_compare,ksize=ksize)
    protein_df = helperfuncs.containments(mat_df=protein_compare,ksize=ksize)
    nt_df.to_csv(f'{WD}/nt_containment{ksize}.csv')
    protein_df.to_csv(f'{WD}/prot_containment{ksize}.csv')

    ### Produce csv file with nt and protein containments with dNdS estimates
    report_df = dnds.report_dNdS(nt_df,protein_df)
    report_df.to_csv(f'{WD}/{output}')

    #create figures of CI analysis

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
