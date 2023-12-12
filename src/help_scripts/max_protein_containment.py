 #!/usr/bin/env python3
"""Return the max C(A,B)"""

import argparse, time, os, subprocess
from dNdS import reportCI
import pandas as pd
pd.set_option('display.max_columns', None)

def main(args):

    ###Arguments
    WD = args.wd
    dnds_constant = args.dnds_constant
    ksize = args.k
    output = args.o

    ### Produce csv file with max nt containment with dNdS estimates
    report_df = reportCI.grab_max_protein_containment_from_containment_csv_file(dnds_constant)
    report_df.to_csv(f'{WD}{output}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Return the max C(A,B)'
    )

    parser.add_argument(
    '--k',
    type=int,
    help = 'Identify the ksize used to produce containment indexes between two sequences.\
    ksize is required to be the same used in both the containment indexes calculated for nucleotide and protein sequences.'
    )

    parser.add_argument(
    '--dnds_constant',
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


'''wd=''
data=pd.read_csv(f'{wd}/test.csv',delimiter=',')
gene_set = set(data['A'].tolist() + data['B'].tolist())
data = data.set_index(['A','B'])
for i in gene_set:
    for j in gene_set:
        if (i,j) in data.index.tolist() and (j,i) in data.index.tolist():
            if data.loc[i].loc[j]['containment_nt'] > data.loc[j].loc[i]['containment_nt']:
                data = data.drop(index=(j, i))
            else:
                data =data.drop(index=(i,j))
        
data.to_csv('max_containment_k.csv')





'''