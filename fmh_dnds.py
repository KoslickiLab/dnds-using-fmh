#!/usr/bin/env python3
"""Approach to estimaating dN/dS ratio of genomes within metagenomic data"""

import argparse
from src import dnds,sourmash_ext
import subprocess
import time

def main(args):
    
    ### PREPARING RUN
    dna_k = args.ksize*3

    if args.mode == "main":
        # Sketch DNA and protein signatures at once using sourmash branchwater plugin manysketch function
        sourmash_ext.run_manysketch(
            fasta_file_csv=args.fasta_input_list,
            ksize=args.ksize,
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

    elif args.mode == "6frametranslate":
         ## Create signatures folder for translate sketches and parallel file
        subprocess.run(f'mkdir {args.directory}/translate_signatures', shell=True, check=True)
        #subprocess.run(f"tail -n +2 {args.fasta_input_list} | grep -v '^$' | cut -d ',' -f 2 | sed 's/^/sourmash sketch translate -f -p k={args.ksize},scaled={args.scaled_input} /' | sed 's/.*\\/\\([^\\/]*\\)\\.fna\\.gz\?/& -o translate_signatures\\/\\1.sig/' > {args.directory}/parallel.txt", shell=True, check=True)
        cmd = (
        f"tail -n +2 {args.fasta_input_list} | grep -v '^$' | cut -d ',' -f 2 "
        f"| sed 's/^/sourmash sketch translate -f -p k={args.ksize},scaled={args.scaled_input} /' "
        r"| sed 's/\([^ ]*\)$/& -o translate_signatures\/\1.sig/' "
        f"> {args.directory}/parallel.txt"
        )
        subprocess.run(cmd, shell=True, check=True)
        subprocess.run(f'parallel < {args.directory}/parallel.txt', shell=True, check=True)
        ## Concat translate sketches
        subprocess.run(f"sourmash signature cat {args.directory}/translate_signatures/ -o {args.directory}/translate.zip", shell=True, check=True)
        ## Run manysketch to produce dna signature
        sourmash_ext.run_manysketch(fasta_file_csv=args.fasta_input_list, ksize=dna_k, scaled=args.scaled_input, cores=args.cores, working_dir=args.directory, molecule='dna')
        ### Run pairwise instead to estimate cfracs
        sourmash_ext.run_pairwise(zipfile=f'{args.directory}/dna.zip',ksize=dna_k,scaled=args.scaled_input,out_csv=f'{args.directory}/results_dna_{dna_k}.csv',molecule='DNA',cores=args.cores, threshold=args.threshold)        
        sourmash_ext.run_pairwise(zipfile=f'{args.directory}/translate.zip',ksize=args.ksize,scaled=args.scaled_input,out_csv=f'{args.directory}/results_translate_{args.ksize}.csv',molecule='protein',cores=args.cores, threshold=args.threshold)
        ### Produce csv file with nt and protein containments with FMH OMEGA estimates
        report_dnds = dnds.report_dNdS_pairwise(f"{args.directory}/results_dna_{dna_k}.csv",f"{args.directory}/results_translate_{args.ksize}.csv",ksize=args.ksize)
        report_dnds.to_csv(f'{args.directory}/fmh_omega_{args.ksize}.csv')

    elif args.mode == "reference":
        # Same sketches, same containment and same dN/dS formula as `main`, for collections of closely
        # related genomes: all pairs are computed from each sketch's differences to one reference sketch.
        # See src/reference_pairwise.py for the identity this relies on.
        from src import reference_pairwise
        if args.reference is None:
            raise ValueError("--mode reference needs --reference, the name of one of the genomes")
        zip_paths = args.sketches
        sketch_start = time.time()
        if not zip_paths and args.sketch_shards > 1:
            zip_paths = reference_pairwise.sketch_in_shards(
                fasta_input_list=args.fasta_input_list,
                ksize=args.ksize,
                scaled=args.scaled_input,
                cores=args.cores,
                directory=args.directory,
                shards=args.sketch_shards
            )
        if not zip_paths:
            sourmash_ext.run_manysketch(
                fasta_file_csv=args.fasta_input_list,
                ksize=args.ksize,
                scaled=args.scaled_input,
                cores=args.cores,
                working_dir=args.directory
            )
            zip_paths = [f'{args.directory}/data.zip']
        sketch_seconds = 0.0 if args.sketches else time.time() - sketch_start
        downsample_to = [int(s) for s in args.downsample_to.split(",")] if args.downsample_to else []
        reference_pairwise.run(
            zip_paths=zip_paths,
            reference_name=args.reference,
            ksize=args.ksize,
            scaled=args.scaled_input,
            directory=args.directory,
            cores=args.cores,
            threshold=args.threshold,
            names_file=args.names,
            downsample_to=downsample_to,
            write_csv=args.write_csv,
            sketch_seconds=sketch_seconds
        )

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

    parser.add_argument(
        '--mode',
        type=str,
        default='main',
        help = 'Enables the use of multithreading from sourmash branchwater plugin and testing translated sequence from sketches. We used --mode 6frametranslate for testing sketch translate, estimations from this mode are currently not reliable for users.\
        Use --mode reference for all pairs of a large collection of closely related genomes: same sketches, containments and dN/dS formula as main, computed from differences to one reference genome (see README).'
    )

    # options of --mode reference (ignored by the other modes)
    parser.add_argument(
        '--reference',
        type=str,
        default=None,
        help = '--mode reference: name (first column of the input csv) of the genome whose sketch the others are compared with. Any genome gives the same numbers; one close to all the others is fastest.'
    )

    parser.add_argument(
        '--sketches',
        type=str,
        nargs='+',
        default=None,
        help = '--mode reference: use these existing sketch zip files (made by this tool or by sourmash scripts manysketch) instead of sketching --fasta_input_list.'
    )

    parser.add_argument(
        '--names',
        type=str,
        default=None,
        help = '--mode reference: file with one genome name per line; only these genomes (and the reference) are compared.'
    )

    parser.add_argument(
        '--downsample_to',
        type=str,
        default=None,
        help = '--mode reference: comma-separated scale factors, coarser than --scaled_input, at which to repeat the comparison by downsampling the same sketches (e.g. 2,5,10).'
    )

    parser.add_argument(
        '--sketch_shards',
        type=int,
        default=1,
        help = '--mode reference: sketch the input csv in this many parts side by side (same sketches, one zip per part). Worth it for thousands of genomes at small scale factors, where writing a single zip is the slow step.'
    )

    parser.add_argument(
        '--write_csv',
        action='store_true',
        help = '--mode reference: also write each result as one fmh_omega_<ksize>.csv, as main mode does (sensible for up to a few million pairs).'
    )


    args = parser.parse_args()

    main(args)
