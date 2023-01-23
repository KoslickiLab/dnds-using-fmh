#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse, time, os, subprocess
from django.urls import path 
from dNdS import findORFs, predictORF, estimatedNdS, reportCI, figures

def main(args):
    """MVC Controller for metagenomic dN/dS estimation
    
    The controller is responsible for:
    -
    -
    """

    genome = args.fasta1
    samples = args.fasta2
    kmers = args.k
    results = args.wd
    
    if args.predict == "orfs":
        """The fasta input file does not have open reading frames identified"""
        findORFs.ORFs_file(genome,results+'ORFs_Fasta1.faa')
        findORFs.ORFs_file(samples,results+'ORFs_Fasta2.faa')
    elif args.predict == "orf":
        """The first fasta file does not have open reading frames identified"""
        print("Getting ORFs for sample input",time.time())
        findORFs.ORFs_file(samples, results+"ORFs_samples.faa")

    if os.path.exists(results+"ORFs_samples.faa"):
        print("Running sourmash sketch...",time.time())
        predictORF.sketch(genome, results+"ORFs_samples.faa",kmer_size=kmers,outputfile1=results+"ref-genome.sig", outputfile2=results+"samples.sig.zip")
        cmd3 = f"unzip {results}samples.sig.zip" #produces SOURMASH-MANIFEST.csv
        subprocess.run(cmd3, stdout=subprocess.PIPE, shell=True)
    else:
        print("File does not exist. Stop analysis.")

    if os.path.exists(results+"SOURMASH-MANIFEST.csv") and os.path.exists(results+"samples.sig.zip") and os.path.exists(results+"ref-genome.sig"):
        print("Running sourmash search...",time.time())
        cmd4 = f"mkdir {results}md5"
        subprocess.run(cmd4, stdout=subprocess.PIPE, shell=True)
        predictORF.search(MANIFEST_CSV_FILE=results+"SOURMASH-MANIFEST.csv", ref_sig=results+"ref-genome.sig", query_sig=results+"samples.sig.zip",output_directory=results+"md5/")
    else:
        print("File does not exist. Stop analysis.")

    if os.path.exists(results+"samples.sig.zip") and os.path.exists(results+"ref-genome.sig"):
        print("Running sourmash prefetch...",time.time())
        predictORF.prefetch(ref_sig=results+"ref-genome.sig", query_sig=results+"samples.sig.zip",kmer_size=kmers)
    else:
        print("File does not exist. Stop analysis.")

def main2(args):

    genome = args.fasta1
    samples = args.fasta2
    kmers = args.k
    results = args.wd

    if args.predict == "orfs":
        """The fasta input file does not have open reading frames identified"""
    #    findORFs.ORFs_file(genome,results+'ORFs_Fasta1.faa')
        findORFs.ORFs_untranslated_file(samples,results+'ORFs_Fasta2.faa')
    elif args.predict == "orf":
        """The first fasta file does not have open reading frames identified"""
        print("Getting ORFs for sample input",time.time())
        findORFs.ORFs_file(samples, results+"ORFs_samples.faa")

    #if os.path.exists(results+"orfs7-70.sig.zip"):
    #    predictORF.prefetch(query_sig=results+"orfs7-70.sig.zip", ref_sig=results+"ref-genome7-70.sig.zip",kmer_size=kmers)
    #else:
    #    print("File does not exist. Stop analysis.")
    #    print(os.path.exists(results+"orfs7-70.sig.zip"))
    #    print(os.path.exists(results+"ref-genome7-70.sig.zip"))

    #reportCI.produce_containment_csv(wd=results)
    #reportCI.prep(data=results+'results_kmers.csv',output=results+"containment.csv")
    #reportCI.analysis_frame01(data=results+"containment.csv",output=results+"CIdict_frame01.pickle")
    #reportCI.analysis_wrongFrames(data=results+"containment.csv",output=results+"CIdict_wrongframes.pickle")

     #create figures of CI analysis
    #figures.CIbox_frame01(data=results+"CIdict_frame01.pickle",wd=results)
    #figures.CIhist(input=results+"containment.csv",wd=results)
    #figures.CIbox_wrongORFs(data=results+"CIdict_wrongframes.pickle",wd=results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'dN/dS estimator for metagenomic data using the containment index between k-mers sets of genomic samples'
    ) 

    parser.add_argument(
        '--fasta1',
        type = str,
        help = 'First fasta file used for dN/dS estimation'
    )

    parser.add_argument(
        '--fasta2',
        type = str,
        help = 'Second fasta file used for dN/dS estimation'
    )

    parser.add_argument(
        '--predict',
        default = 'orf',
        choices=['', 'orf', 'orfs'],
        help = 'Flagged when fasta files need open reading frames prediction for all fasta files (orfs), for the first fasta file (orfs1), or for the second fasta file (orfs2)'
    )

    parser.add_argument(
        '--k',
        default = [7],
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
        default = 'dNdS_estimates.txt',
        help = 'Output file name for dN/dS report'
    )

    parser.add_argument(
        '--wd',
        help = 'Output files from sourmash program'
    )

    args = parser.parse_args()

    main2(args)
