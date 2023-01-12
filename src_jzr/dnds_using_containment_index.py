#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse

from dNdS import findORFs, predictORF, estimatedNdS

def main(args):
    """MVC Controller for metagenomic dN/dS estimation
    
    The controller is responsible for:
    -
    -
    """

    fasta1 = args.fasta1
    fasta2 = args.fasta2

    if args.predict == "orfs":
        """The fasta input file does not have open reading frames identified"""
        translate_file(fasta1_file,'ORFs_Fasta1.faa')
        translate_file(fasta2_file,'ORFs_Fasta2.faa')

    elif args.predict == "orfs1":
        """The first fasta file does not have open reading frames identified"""
        translate_file(fasta1_file, 'ORFs_Fasta1.faa')
    
    elif args.predict == "orfs2":
        """The second fasta file does not have open reading frames identified"""
        translate_file(fasta2_file, 'ORFs_Fasta2.faa')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'dN/dS estimator for metagenomic data using the containment index between k-mers sets of genomic samples'
    ) 

    parser.add_argument(
        '--fasta1',
        nargs = 1,
        help = 'First fasta file used for dN/dS estimation'
    )

    parser.add_argument(
        '--fasta2',
        nargs = 1,
        help = 'Second fasta file used for dN/dS estimation'
    )

    parser.add_argument(
        '--predict',
        default = '',
        choices=['', 'orfs', 'orfs1', 'orfs2'],
        help = 'Flagged when fasta files need open reading frames prediction for all fasta files (orfs), for the first fasta file (orfs1), or for the second fasta file (orfs2)'
    )

    parser.add_argument(
        '--k',
        default = 7,
        type = int,
        help = 'K-mer size for containment index, a parameter used in Sourmash sketch.'
    )

    parser.add_argument(
        '--scaled1',
        default = 100,
        type = int,
        help = 'Scale first fasta input file, a parameter used in Sourmash sketch'
    )

    parser.add_argument(
        '--scaled2',
        default = 1,
        type = int,
        help = 'Scale second fasta input file, a parameter used in Sourmash sketch'
    )

    parser.add_argument(
        '--Tbp',
        default = 1,
        type = int,
        help = 'threshold-bp, a parameter used in Sourmash prefetch'
    )

    parser.add_argument(
        '--output',
        default = 'XYZ',
        help = 'Output file name for dN/dS report'
    )

    args = parser.parse_args()

    main(args)
