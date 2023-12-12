 #!/usr/bin/env python3
"""Return the max C(A,B)"""

import argparse, time, os, subprocess
from dNdS import reportCI
import pandas as pd
pd.set_option('display.max_columns', None)

def main(args):

    ###Arguments
    WD = args.wd
    max_dnds_constant = args.max_dnds_constant
    output = args.o

    ### Produce csv file with max nt containment with dNdS estimates
    report_df = reportCI.change_column_names(max_dnds_constant)
    report_df.to_csv(f'{WD}{output}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Return the max C(A,B)'
    )

    parser.add_argument(
    '--max_dnds_constant',
    help = 'Input of dnds_constant_k.csv filefile.\
    This file is a pairwise matrix produced from sourmash compare that includes containment indexes between nucleotide sequences.'
    )

    parser.add_argument(
    '--o',
    help = 'Output CSV file with max conainment nt and dN/dS ratio estimates between sequences evaluated. The CSV file reports sequence A, sequence B,\
    ksize, and containment index'
    )

    parser.add_argument(
    '--wd',
    help = 'Output directory for CSV file with dN/dS ratio estimates between sequences evaluated.'
    )

    args = parser.parse_args()

    main(args)

