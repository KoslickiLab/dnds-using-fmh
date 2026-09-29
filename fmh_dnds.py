#!/usr/bin/env python3
"""Approach to estimaating dN/dS ratio of genomes within metagenomic data"""

import argparse
from src import dnds,sourmash_ext
import subprocess

def main(args):
    
    ### PREPARING RUN
    dna_k = args.ksize*3

    # Sketch DNA and protein signatures at once using sourmash branchwater plugin manysketch function
    sourmash_ext.run_manysketch(
        fasta_file_csv=args.fasta_input_list,
        kaa=args.ksize,
        scaled=args.scaled_input,
        cores=args.cores,
        working_dir=args.directory
    )

    # Estimate containments between sketches before fmh dnds estimations
    sourmash_ext.run_pairwise(
        zipfile=f'{args.directory}/data.zip',
        k = dna_k,
        scaled=args.scaled_input,
        out_csv=f'{args.directory}/results_dna_{dna_k}.csv',
        molecule='DNA',
        cores=args.cores,
        threshold=args.threshold
    )
    sourmash_ext.run_pairwise(
        zipfile=f'{args.directory}/data.zip',
        k=args.ksize,
        scaled=args.scaled_input,
        out_csv=f'{args.directory}/results_protein_{args.ksize}.csv',
        molecule='protein',
        cores=args.cores,
        threshold=args.threshold
    )

    # Estimate dN/dS and report
    report_dnds = dnds.report_dNdS_pairwise(
        f"{args.directory}/results_dna_{dna_k}.csv",
        f"{args.directory}/results_protein_{args.ksize}.csv",
        ksize=args.ksize
    )
    report_dnds.to_csv(f'{args.directory}/fmh_omega_{args.ksize}.csv')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'dN/dS estimator for metagenomic data using the containment index between k-mer sets of genomic samples'
    )

    parser.add_argument(
        '--fasta_input_list',
        nargs='?',
        const='arg_was_not_given',
        help = 'Input csv file that contains fasta files for sketching.\
        This csv file follows sourmash scripts example, where in the first column is name, second column is dna fasta filename, and third column is protein fasta filename.'
    )


    parser.add_argument(
        '--ksize',
        type=int,
        default=7,
        help = 'Identify a ksize used to produce sketches. Specifically, here it refers to the protein ksize.\
        ksize is required to be the same used in both the containment indexes calculated for nucleotide and protein sequences.'
    )

    parser.add_argument(
        '--scaled_input',
        type=int,
        default=500,
        help = 'Identify a scaled factor for signature sketches.\
        Use a scale factor of at least 10 for thousands of genomes. The scale factor represents the compression level used (higher=more compression).'
    )

    parser.add_argument(
        '--directory',
        type=str,
        help = 'Output directory for FMH Omega estimation.'
    )

    parser.add_argument(
        '--cores',
        type=int,
        default=1,
        help = 'Set total cores. Use anything above 100 when using thousands of genomes.'
    )

    parser.add_argument(
        '--threshold',
        type=float,
        default=0.05,
        nargs='?',
        const='arg_was_not_given',
        help = 'Set containment threshold for sourmash plugin branchwater commands. In short, dN/dS values will not be calculated if the containment index is below this threshold (i.e. too distant of genomic sequences; not similar enough to compare).'
    )    

    args = parser.parse_args()

    main(args)
