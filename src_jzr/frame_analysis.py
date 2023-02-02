#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse, time, os, subprocess
#from django.urls import path 
#from dNdS import findORFs, predictORF, estimatedNdS, reportCI
from dNdS import reportCI,figures

def main(args):
    """MVC Controller for metagenomic dN/dS estimation
    
    The controller is responsible for:
    -
    -
    """

    results = args.wd

    ### The following functions are to evaluate and produce figures for frame prediction analysis
    # report highest and second highest containment index in dictionary pickle file
    #reportCI.produce_containment_csv(wd=results)
    #reportCI.prep(data=results+'results_kmers.csv',output=results+"containment.csv")
    #reportCI.analysis_frame1(data=results+"containment.csv",output=results+"CIdict.pickle")
    reportCI.extract_frame1(data=results+"containment.csv",output=results+"frame1_CI.csv")
    reportCI.extract_frameX(data=results+"containment.csv",output=results+"framex_CI.csv")

    #create figures of CI analysis
    figures.CIbox_frames(frame_1data=results+"frame1_CI.csv",frame_xdata=results+"framex_CI.csv",wd=results)
    #figures.CIhist(data=results+"containment.csv",wd=results)
    #figures.CIboxplots_frames(data=results+"containment.csv",wd=results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'dN/dS estimator for metagenomic data using the containment index between k-mers sets of genomic samples'
    ) 

    #parser.add_argument(
    #    '--ref_fasta',
    #    type = str,
    #    help = 'First fasta file used for dN/dS estimation'
    #)

    #parser.add_argument(
    #    '--query_fasta',
    #    type = str,
    #    help = 'Second fasta file used for dN/dS estimation'
    #)

    #parser.add_argument(
    #    '--predict',
    #    default = 'orf',
    #    choices=['', 'orf', 'orfs','frame'],
    #    help = 'Flagged when fasta files need open reading frames prediction for all fasta files (orfs), for the first fasta file (orfs1), or for the second fasta file (orfs2)'
    #)

    #parser.add_argument(
    #    '--k',
    #    default = [7],
    #    help = 'K-mer size for containment index, a parameter used in Sourmash sketch.'
    #)

    #parser.add_argument(
    #    '--scaled1',
    #    default = 100,
    #    type = int,
    #    help = 'Scale first fasta input file, a parameter used in Sourmash sketch'
    #)

    #parser.add_argument(
    #    '--scaled2',
    #    default = 1,
    #    type = int,
    #    help = 'Scale second fasta input file, a parameter used in Sourmash sketch'
    #)

    #parser.add_argument(
    #    '--Tbp',
    #    default = 1,
    #    type = int,
    #    help = 'threshold-bp, a parameter used in Sourmash prefetch'
    #)

    #parser.add_argument(
    #    '--output',
    #    default = 'dNdS_estimates.txt',
    #    help = 'Output file name for dN/dS report'
    #)

    parser.add_argument(
        '--wd',
        help = 'Output files from sourmash program'
    )

#    parser.add_argument(
#        '--moltype',
#        default = 'protein',
#        choices = ['protein','dna'],
#        help = 'Indicate if protein or DNA sequences are being used'
#    )

#    parser.add_argument(
#        '--translate',
#        default = 'no',
#        choices = ['yes','no'],
#        help = 'Indicate we are translating both ref and query sequences.'
#    )

    #parser.add_argument(
    #    '--analyze',
    #    default = 'frame_predict',
    #    choices = ['frame_predict','dna','protein'],
    #    help = 'type of sourmash analysis'
    #)

    args = parser.parse_args()

    main(args)
