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

    genome = args.ref_fasta #reference
    samples = args.query_fasta #queries
    kmers = args.k
    results = args.wd
    sd1 = args.scaled1
    analysis_type = args.analyze
#    mltyp = args.moltype
#    trnslt = args.translate
    
    #if args.predict == "orfs":
    #    """The fasta input file does not have open reading frames identified"""
    #    findORFs.ORFs_file(genome,results+'ORFs_ref_fasta.faa')
    #    findORFs.ORFs_file(samples,results+'ORFs_query_fasta.faa')
    #elif args.predict == "orf":
    #    """The first fasta file does not have open reading frames identified"""
    #    print("Getting ORFs for sample input",time.time())
    #    findORFs.ORFs_file(samples, results+"ORFs_samples.faa")
#    if args.moltype == 'protein':
#        print('Analyzing protein sequences')
    if args.analyze == 'frame_predict':
            """Create fasta file with six reading frames of each sequence"""
            print("Obtain six reading frames for each query")
            findORFs.reading_frames_file(samples, results+"query_frames.faa")
            query = results+"query_frames.faa"
#    elif args.moltype == 'dna':
    elif args.analyze == "dna" or args.analyze == "protein":
            print('Analyzing DNA sequences')
            query = samples

    print('Ready to sketch!')
    if os.path.exists(query):
        print("Running sourmash sketch...",time.time())
#        predictORF.sketch(genome, query,kmer_size=kmers,ref_output=results+"ref-genome.sig", query_output=results+"queries.sig.zip",moltype=mltyp,translate=trnslt)
        predictORF.sketch(genome, query,kmer_size=kmers,ref_output=results+"ref-genome.sig", query_output=results+"queries.sig.zip",analysis=analysis_type)
        #cmd3 = f"unzip {results}queries.sig.zip" #produces SOURMASH-MANIFEST.csv How do I output to a directory
        #subprocess.run(cmd3, stdout=subprocess.PIPE, shell=True)
    #if os.path.exists(results+samples):
    #    print("Running sourmash sketch...",time.time())
    #    predictORF.sketch(genome, samples, kmer_size=kmers, scaledfile1=sd1, ref_output=results+"ref-genome.sig", query_output=results+"queries.sig.zip")
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
#        predictORF.prefetch(query_sig="queries.sig.zip", ref_sig="ref-genome.sig",kmer_size=kmers,wd=results,moltype=mltyp)
        predictORF.prefetch(query_sig="queries.sig.zip", ref_sig="ref-genome.sig",kmer_size=kmers,wd=results,analysis=analysis_type)
    else:
        print("File does not exist. Stop analysis.")
        print(os.path.exists(results+"queries.sig.zip"))
        print(os.path.exists(results+"ref-genome.sig"))

    ### The following functions are to evaluate and produce figures for frame prediction analysis
    # report highest and second highest containment index in dictionary pickle file
    #reportCI.produce_containment_csv(wd=results)
    #reportCI.prep(data=results+'results_kmers.csv',output=results+"containment.csv")
    #reportCI.analysis_frame01(data=results+"containment.csv",output=results+"CIdict.pickle")

    #create figures of CI analysis
    #figures.CIhist(data=results+"CIdict.pickle",wd=results)
    #figures.CIbox_frame01(input=results+"containment.csv",wd=results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'dN/dS estimator for metagenomic data using the containment index between k-mers sets of genomic samples'
    ) 

    parser.add_argument(
        '--ref_fasta',
        type = str,
        help = 'First fasta file used for dN/dS estimation'
    )

    parser.add_argument(
        '--query_fasta',
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

    parser.add_argument(
        '--analyze',
        default = 'frame_predict',
        choices = ['frame_predict','dna','protein'],
        help = 'type of sourmash analysis'
    )

    args = parser.parse_args()

    main(args)
