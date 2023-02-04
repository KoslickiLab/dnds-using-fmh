#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse, time, os, subprocess
#from django.urls import path 
#from dNdS import findORFs, predictORF, estimatedNdS, reportCI
from dNdS import reportCI,reportdNdS

def main(args):
    """MVC Controller for metagenomic dN/dS estimation"""

    nt_compare = args.nt
    protein_compare = args.protein
    ksize = args.k
    output = args.o

    ### Obtain containment from matrix files produced by sourmash compare
    print(nt_compare)
    nt_df = reportCI.grab_containment_from_mat(mat_df=nt_compare,ksize=ksize)
    protein_df = reportCI.grab_containment_from_mat(mat_df=protein_compare,ksize=ksize)

    ### Produce csv file with nt and protein containments with dNdS estimates
    report_df = reportdNdS.report_dNdS(nt_df,protein_df)
    report_df.to_csv(output)

    #create figures of CI analysis

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'dN/dS estimator for metagenomic data using the containment index between k-mers sets of genomic samples'
    )

    parser.add_argument(
        '--k',
        type=int,
        help = 'Output files from sourmash program'
    )

    parser.add_argument(
        '--nt',
        help = 'input nt csv file'
    )

    parser.add_argument(
        '--protein',
        help = 'Input protein csv file'
    )

    parser.add_argument(
    '--o',
    help = 'Output csv file with dNdS estimates'
    )

    args = parser.parse_args()

    main(args)
