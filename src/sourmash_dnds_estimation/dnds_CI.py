#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse, time, os, subprocess
from dNdS import reportCI,CfracdNdS

def main(args):

    ###Arguments
    nt_compare = args.nt
    protein_compare = args.protein
    ksize = args.k
    output = args.o

    ### Obtain containment from matrix files produced by sourmash compare
    nt_df = reportCI.grab_containment_from_mat_ground_truth(mat_df=nt_compare,ksize=ksize)
    protein_df = reportCI.grab_containment_from_mat_ground_truth(mat_df=protein_compare,ksize=ksize)

    ### Produce csv file with nt and protein containments with dNdS estimates
    report_df = CfracdNdS.report_dNdS(nt_df,protein_df)
    report_df.to_csv(output)

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

    args = parser.parse_args()

    main(args)
