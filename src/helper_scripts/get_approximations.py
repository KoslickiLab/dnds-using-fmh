#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse
from fmh_omega import dnds
import pandas as pd
import subprocess
import time
import multiprocessing as mp
import numpy as np
import math

def main(args):
    
    ### ARGUMENTS
    data = args.input
    k = args.ksize
    on = args.outname

    df = pd.read_csv(f'{data}').rename(columns={'DNA_Cfrac':'DNA_max_Cfrac','AA_Cfrac':'AA_max_Cfrac'})
    df['ANI_approx'] = dnds.ANI_approx(df['DNA_max_Cfrac'],k) #this function gets dna ksize for you
    df['AAI_approx'] = dnds.AAI_approx(df['AA_max_Cfrac'],k)

    df.to_csv(f'{on}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'output information with ANI and AAI approximations'
    )

    parser.add_argument(
        '--input',
        help = 'Input txt file that contains fasta files for sketching.\
        In first column, have fasta file name and second column have nickname for each fasta files.'
    )

    parser.add_argument(
        '--ksize',
        type=int,
        help = 'Identify a ksize used to produce sketches.\
        ksize is required to be the same used in both the containment indexes calculated for nucleotide and protein sequences.'
    )

    parser.add_argument(
        '--outname',
        help = 'Name your study for FMH OMEGA Estimations. Prefix for output files.'
    )

    args = parser.parse_args()

    main(args)

