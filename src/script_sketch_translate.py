#!/usr/bin/env python3
"""This script creates a parallel job and runs it to sketch translate on thousands of genomes"""

import argparse
import subprocess

def main(args):
    
    ### get dna ksize
    dna_k = args.ksize*3

    ## Create signatures folder for translate sketches and parallel file
    subprocess.run(f'mkdir {args.directory}/translate_signatures', shell=True, check=True)
    subprocess.run(f"tail -n +2 {args.fasta_input_list} | cut -d ',' -f 2 | sed 's/^/sourmash sketch translate -f -p k={args.ksize},scaled={args.scaled_input} /' | sed 's/.*\\/\\([^\\/]*\\)\\.fna\\.gz/& -o translate_signatures\\/\\1.sig.gz/' > {args.directory}/parallel.txt", shell=True, check=True)
    subprocess.run(f'parallel < {args.directory}/parallel.txt', shell=True, check=True)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'This script serves to run parallel for sketch translate on thousands of genomes'
    )

    parser.add_argument(
        '--fasta_input_list',
        nargs='?',
        const='arg_was_not_given',
        help = 'Input csv file that contains fasta files for sketching.\
        This csv file follows sourmash scripts example, where in the first column is name, second column is dna fasta filename, and third column is protein fasta filename.'
    )

    parser.add_argument(
        '--dna_fasta',
        nargs='?',
        const='arg_was_not_given',
        help = 'Input filename that contains dna entries for sketching using --singleton parameter.'
    )

    parser.add_argument(
        '--protein_fasta',
        nargs='?',
        const='arg_was_not_given',
        help = 'Input filename that contains protein entries for sketching using --singleton parameter.'
    )

    parser.add_argument(
        '--ksize',
        type=int,
        help = 'Identify a ksize used to produce sketches. Specifically, here it refers to the protein ksize.\
        ksize is required to be the same used in both the containment indexes calculated for nucleotide and protein sequences.'
    )

    parser.add_argument(
        '--scaled_input',
        type=int,
        help = 'Identify a scaled factor for signature sketches.\
        E.g., A scaled factor = 1 will include all k-mers in final signature.'
    )

    parser.add_argument(
        '--mode',
        type=str,
        help = 'Identify mode to run fmh_omega as sngl, mult, bwmult, bwpair'
    )

    parser.add_argument(
        '--translate',
        type=str,
        help='indicate yes or no to translate coding sequence'
    )

    parser.add_argument(
        '--outname',
        help = 'Name your study for FMH OMEGA Estimations. Prefix for output files.'
    )

    parser.add_argument(
        '--directory',
        type=str,
        help = 'Output directory for FMH Omega estimation.'
    )

    parser.add_argument(
        '--cores',
        type=int,
        help = 'Set total cores. Use anything above 100 when using thousands of genomes.'
    )

    parser.add_argument(
        '--threshold',
        type=float,
        default=0.0,
        nargs='?',
        const='arg_was_not_given',
        help = 'Set containment threshold for sourmash plugin branchwater commands.'
    )    

    args = parser.parse_args()

    main(args)

