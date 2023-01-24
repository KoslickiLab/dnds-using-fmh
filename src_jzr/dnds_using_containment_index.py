#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse, time, os, subprocess
from django.urls import path 
from dNdS import findORFs, predictORF, estimatedNdS, reportCI

def main(args):
    """MVC Controller for metagenomic dN/dS estimation
    
    The controller is responsible for:
    -
    -
    """

    genome = args.fasta1
    samples = args.fasta2 #queries
    kmers = args.k
    results = args.wd
    sd1 = args.scaled1
    
    if args.predict == "orfs":
        """The fasta input file does not have open reading frames identified"""
        findORFs.ORFs_file(genome,results+'ORFs_Fasta1.faa')
        findORFs.ORFs_file(samples,results+'ORFs_Fasta2.faa')
    elif args.predict == "orf":
        """The first fasta file does not have open reading frames identified"""
        print("Getting ORFs for sample input",time.time())
        findORFs.ORFs_file(samples, results+"ORFs_samples.faa")
    elif args.predict == "frame":
        """Create fasta file with six reading frames of each sequence"""
        print("Obtain six reading frames for each query")
        findORFs.reading_frames_file(samples, results+"query_frames.faa")
        query = results+"query_frames.faa"

    print('Ready to sketch!')
    if os.path.exists(query):
        print("Running sourmash sketch...",time.time())
        predictORF.sketch(genome, query,kmer_size=kmers,outputfile1=results+"ref-genome.sig", outputfile2=results+"queries.sig.zip")
        cmd3 = f"unzip {results}queries.sig.zip" #produces SOURMASH-MANIFEST.csv
        subprocess.run(cmd3, stdout=subprocess.PIPE, shell=True)
    #if os.path.exists(results+samples):
    #    print("Running sourmash sketch...",time.time())
    #    predictORF.sketch(genome, samples, kmer_size=kmers, scaledfile1=sd1, outputfile1=results+"ref-genome.sig", outputfile2=results+"queries.sig.zip")
    #    cmd3 = f"unzip {results}samples.sig.zip" #produces SOURMASH-MANIFEST.csv
    #    subprocess.run(cmd3, stdout=subprocess.PIPE, shell=True)
    else:
        print("File does not exist. Stop analysis.")

    #this is not part of pipeline
    #if os.path.exists(results+"SOURMASH-MANIFEST.csv") and os.path.exists(results+"samples.sig.zip") and os.path.exists(results+"ref-genome.sig"):
    #    print("Running sourmash search...",time.time())
    #    cmd4 = f"mkdir {results}md5"
    #    subprocess.run(cmd4, stdout=subprocess.PIPE, shell=True)
    #    predictORF.search(MANIFEST_CSV_FILE=results+"SOURMASH-MANIFEST.csv", ref_sig=results+"ref-genome.sig", query_sig=results+"samples.sig.zip",output_directory=results+"md5/")
    #else:
    #    print("File does not exist. Stop analysis.")

    print('Ready to prefetch!')
    if os.path.exists(results+"queries.sig.zip"):
        print("Running sourmash prefetch...",time.time())
        predictORF.prefetch(query_sig="queries.sig.zip", ref_sig="ref-genome.sig",kmer_size=kmers,wd=results)
    else:
        print("File does not exist. Stop analysis.")
        print(os.path.exists(results+"queries.sig.zip"))
        print(os.path.exists(results+"ref-genome.sig"))

    # report highest and second highest containment index in dictionary pickle file
    reportCI.produce_containment_csv(wd=results)
    reportCI.prep(data=results+'results_kmers.csv',output=results+"containment.csv")
    reportCI.analysis_frame01(data=results+"containment.csv",output=results+"CIdict.pickle")

    #create figures of CI analysis
    figures.CIhist(data=results+"CIdict.pickle",wd=results)
    figures.CIbox_frame01(input=results+"containment.csv",wd=results)

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
        choices=['', 'orf', 'orfs','frame'],
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

    main(args)
