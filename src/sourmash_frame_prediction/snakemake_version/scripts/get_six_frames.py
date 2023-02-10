#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse, time, os, subprocess
from frame_predict import findORfs

def main(args):
    """MVC Controller for metagenomic dN/dS estimation
    
    The controller is responsible for:
    -
    -
    """

    results = args.wd

    """Create fasta file with six reading frames of each sequence"""
    print("Obtain six reading frames for each query")
    findORFs.reading_frames_file(samples, results+"query_frames.faa")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'dN/dS estimator for metagenomic data using the containment index between k-mers sets of genomic samples'
    ) 

    parser.add_argument(
        '--query_fasta',
        type = str,
        help = 'Second fasta file used for dN/dS estimation'
    )

    parser.add_argument(
        '--output',
        default = 'dNdS_estimates.txt',
        help = 'Output file name for dN/dS report'
    )

    parser.add_argument(
        '--wd',
        help = 'Output files from sourmash program'
    )

    args = parser.parse_args()

    main(args)
