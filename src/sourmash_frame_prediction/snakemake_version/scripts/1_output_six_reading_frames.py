#!/usr/bin/env python3
"""Obtain six reading frames from nucleotide protein coding sequences"""

import argparse, time, os, subprocess
#from django.urls import path 
#from dNdS import findORFs, predictORF, estimatedNdS, reportCI
from frame_predict import findORFs

def main(args):
    """arguments"""
    fasta_file = args.nucleotide_fasta
    output_file = args.output

    """Create frame .faa file"""
    findORFs.reading_frames_file(INFILE=fasta_file,OUTPUT_FILENAME=output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Output a new fasta file with six reading frames of each sequence, where the first reading frame is the correct one'
    ) 

    parser.add_argument(
        '--nucleotide_fasta',
        type = str,
        help = 'Input a fasta file that includes nucleotide fasta sequences of protein coding sequences'
    )

    parser.add_argument(
        '--output',
        default = 'dNdS_estimates.txt',
        help = 'Output file name for dN/dS report'
    )

    args = parser.parse_args()

    main(args)
