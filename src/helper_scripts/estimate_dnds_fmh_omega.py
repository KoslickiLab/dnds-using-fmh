#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data. THis script only estimates dn/dS from FMH. You need containment estimates from sourmash to run this"""

import argparse
from fmh_omega import dnds

def main(args):
    ### Produce csv file with nt and protein containments with FMH OMEGA estimates
    report_dnds = dnds.report_dNdS_pairwise(f"{args.dna_containment}",f"{args.protein_containment}",ksize=args.ksize)
    report_dnds.to_csv(f'{args.directory}/fmh_omega_{args.ksize}.csv')
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'dN/dS estimator for metagenomic data using the containment index between k-mer sets of genomic samples'
    )

    parser.add_argument(
        '--dna_containment',
        type=str,
        help = 'Input filename that contains dna containment from sourmash pariwise.'
    )

    parser.add_argument(
        '--protein_containment',
        type=str,
        help = 'Input filename that contains protein containment from sourmash pariwise.'
    )

    parser.add_argument(
        '--ksize',
        type=int,
        help = 'Identify a ksize used when producing sketches using sourmash. Specifically, here it refers to the protein ksize.\
        ksize is required to be the same used in both the containment indexes calculated for nucleotide and protein sequences.'
    )

    parser.add_argument(
        '--directory',
        type=str,
        help = 'Output directory for FMH Omega estimation.'
    )

    args = parser.parse_args()

    main(args)

