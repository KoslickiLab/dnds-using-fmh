#!/usr/bin/env python3
"""Obtain six reading frames from nucleotide protein coding sequences"""

import argparse, time, os, subprocess
#from django.urls import path 
#from dNdS import findORFs, predictORF, estimatedNdS, reportCI
from frame_predict import findORfs,reportCI,figures

def main(args):
    """arguments"""
    fasta_file = args.nucleotide_fasta
    output_file = args.output

    """Create frame .faa file"""
    findORFs.reading_frames_file(INFILE=fasta_file,OUTPUT_FILENAME=output_file)

    ### The following functions are to evaluate and produce figures for frame prediction analysis
    
    # report highest and second highest containment index in dictionary pickle file
    #reportCI.produce_containment_csv(wd=results)
    #reportCI.prep(data=results+'results_kmers.csv',output=results+"containment.csv")
    #reportCI.analysis_frame1(data=results+"containment.csv",output=results+"CIdict.pickle")
    
    #### report CI for frames 1 and all other frames
    #reportCI.extract_frame1(data=results+"containment.csv",output=results+"frame1_CI.csv")
    #reportCI.extract_frameX(data=results+"containment.csv",output=results+"framex_CI.csv")

    #create figures of CI analysis
    #figures.CIbox_frames(frame_1data=results+"frame1_CI.csv",frame_xdata=results+"framex_CI.csv",wd=results) #all other frames not 1
    #figures.CIbox_frames(frame_1data=results+"frame1_CI.csv",frame_xdata=results+"second_highest_containment.csv",wd=results,output='boxplot_frame1_vs_second_highest.jpeg') #second highest
    #figures.CIhist(data=results+"containment.csv",wd=results)
    #figures.CIboxplots_frames(data=results+"containment.csv",wd=results)

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
