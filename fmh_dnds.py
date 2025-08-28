#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse
from fmh_dnds import helperfuncs,dnds,sourmash_ext
import subprocess
#import time
#import multiprocessing as mp
#import numpy as np
#import math

def main(args):
    
    ### ARGUMENTS
    #dna_fasta = args.fasta_input_list
    #k = args.ksize
    #s = args.scaled_input
    #on = args.outname
    #wd = args.directory
    #m =args.mode
    #translate_cds=args.translate
    #total_cores =args.cores
    #thresh =args.threshold

    ### PREPARING RUN
    dna_k = args.ksize*3


    ### RUN WHEN NOT USING SOURMASH BRANCHWATER PLUGIN
    #if m == "bwmult" and m !="bwpair":
    if args.mode == "mult":
        #store file information in lists
        fastn_files=[] #dna fasta file
        fasta_files=[] #protein fasta files
        for files in [line.strip() for line in open(f'{args.fasta_input_list}', 'r')]:
            if files != '' and 'genome_filename' not in files:
                fastn_filename = files.split(',')[1]
                name = files.split(',')[1].split('.')[0]
                fastn_files.append(fastn_filename)
                fasta_files.append(f'{name}.translated.fasta')
        if args.translate == 'yes':
        #Translate before sketching protein
            for pos in range(len(fastn_files)):
                helperfuncs.translate_CDS(cds_fasta=f'{fastn_files[pos]}', out_name=f'{fasta_files[pos]}')

        #Create signature directory
        subprocess.run(f'mkdir {args.directory}/signatures', shell=True, check=True)
        #Sketch signatures
        #fastn_lst = f'{args.directory}/'+f' {args.directory}/'.join(fastn_files)
        #fasta_lst = f' '.join(fasta_files)
        
        dna_sig_outname = f'{args.directory}/signatures/{args.outname}.dna.sig.gzip'
        prot_sig_outname = f'{args.directory}/signatures/{args.outname}.protein.sig.gzip'

        if m == 'mult': #input are fastn files
            sourmash_ext.sketch_genome_dna(fasta=args.fasta_input_list, ksize=dna_k, scaled=args.scaled_input, out_sigfile=dna_sig_outname, multiple={m})
            sourmash_ext.sketch_genome_protein(fasta=args.fasta_input_list, ksize=args.ksize, scaled=args.scaled_input, out_sigfile=prot_sig_outname, multiple={m})
        elif m == 'sngl':
            sourmash_ext.sketch_genome_dna(fasta=args.fasta_input_list, ksize=dna_k, scaled=args.scaled_input, out_sigfile=dna_sig_outname)
            sourmash_ext.sketch_genome_protein(fasta=args.fasta_input_list, ksize=args.ksize, scaled=args.scaled_input, out_sigfile=prot_sig_outname)
        # Create compare directory
        subprocess.run(f'mkdir {args.directory}/compare_dna', shell=True, check=True)
        subprocess.run(f'mkdir {args.directory}/compare_protein', shell=True, check=True)
        ### run sourmash comapre for cfracs
        sourmash_ext.compare_signatures(ref=f'{args.directory}/signatures/{args.outname}.dna.sig.gzip', query=f'{args.directory}/signatures/{args.outname}.dna.sig.gzip', ksize=dna_k, molecule='dna', working_dir=f'{args.directory}')
        sourmash_ext.compare_signatures(ref=f'{args.directory}/signatures/{args.outname}.protein.sig.gzip', query=f'{args.directory}/signatures/{args.outname}.protein.sig.gzip', ksize=args.ksize, molecule='protein', working_dir=f'{args.directory}')
        #Report containments
        nt_df = helperfuncs.containments(mat_df=f'{args.directory}/compare.dna.{dna_k}.csv',ksize=int(args.ksize),multiple=m)
        protein_df = helperfuncs.containments(mat_df=f'{args.directory}/compare.protein.{args.ksize}.csv',ksize=int(args.ksize),multiple=m)
        ### Produce csv file with nt and protein containments with FMH OMEGA estimates
        report_df = dnds.report_dNdS(nt_df,protein_df)
        report_df.to_csv(f'{args.directory}/fmh_omega_{k}.csv')

    ### RUN WHEN USING SOURMASH BRANCHWATER PLUGIN
    # Run when we are using raw reads or entire genomes
    elif args.mode == "translate" or args.mode == "sngl_translate":
        if args.mode == "translate":
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

        elif args.mode == "sngl_translate":
            sourmash_ext.sketch_genome_dna(fasta=args.dna_fasta, ksize=dna_k, scaled=args.scaled_input, out_sigfile=f'{args.directory}/dna.zip')
            sourmash_ext.sketch_genome_translate(fasta=args.dna_fasta, ksize=args.ksize, scaled=args.scaled_input, out_sigfile=f'{args.directory}/translate.zip')
            ### Run pairwise instead to estimate cfracs
            sourmash_ext.compare_signatures(ref=f'{args.directory}/dna.zip',query=f'{args.directory}/dna.zip',ksize=dna_k,molecule='dna',working_dir=args.directory)
            sourmash_ext.compare_signatures(ref=f'{args.directory}/translate.zip',query=f'{args.directory}/translate.zip',ksize=args.ksize,molecule='protein',working_dir=args.directory)
            #Report containments
            nt_df = helperfuncs.extract_containment_matrix(mat_csv=f'{args.directory}/compare.dna.{dna_k}.csv')
            protein_df = helperfuncs.extract_containment_matrix(mat_csv=f'{args.directory}/compare.protein.{args.ksize}.csv')
            ### Produce csv file with nt and protein containments with FMH OMEGA estimates
            report_df = dnds.report_dNdS_6frame(nt_df,protein_df,ksize=args.ksize) #constant is incorrected
            report_df.to_csv(f'{args.directory}/fmh_omega_{args.ksize}.csv') #does not include p_nt_mut nor p_no_mut
            #report_df.to_csv(f'{args.directory}/fmh_omega_{args.ksize}_with_additional_calculations.csv')
        
    elif args.mode == "bwmult" or args.mode == "bwpair" or args.mode == "sngl":
        ###get total expected signatures
        total_num_signatures=-1
        with open(f'{args.fasta_input_list}') as infp:
            for line in infp:
                if line.strip():
                    total_num_signatures += 1
        sourmash_ext.run_manysketch(fasta_file_csv=args.fasta_input_list, ksize=args.ksize, scaled=args.scaled_input, cores=args.cores, working_dir=args.directory)
        if args.mode == "bwmult":
            ### Run multisearch to estimate cfracs
            sourmash_ext.run_multisearch(ref_zipfile=f'{args.directory}/dna.zip',query_zipfile=f'{args.directory}/dna.zip',ksize=dna_k,scaled=args.scaled_input,out_csv=f'{args.directory}/results_dna_{dna_k}.csv',molecule='DNA',cores=args.cores)
            sourmash_ext.run_multisearch(ref_zipfile=f'{args.directory}/protein.zip',query_zipfile=f'{args.directory}/protein.zip',ksize=args.ksize,scaled=args.scaled_input,out_csv=f'{args.directory}/results_protein_{args.ksize}.csv',molecule='protein',cores=args.cores)
            ### Produce csv file with nt and protein containments with FMH OMEGA estimates
            report_dnds = dnds.report_dNdS_multisearch(f"{args.directory}/results_dna_{dna_k}.csv",f"{args.directory}/results_protein_{args.ksize}.csv",ksize=args.ksize)
            report_dnds.to_csv(f'{args.directory}/fmh_omega_{args.ksize}.csv')
        elif args.mode == "bwpair":
            ### Run pairwise instead to estimate cfracs
            sourmash_ext.run_pairwise(zipfile=f'{args.directory}/data.zip',ksize=dna_k,scaled=args.scaled_input,out_csv=f'{args.directory}/results_dna_{dna_k}.csv',molecule='DNA',cores=args.cores, threshold=args.threshold)
            sourmash_ext.run_pairwise(zipfile=f'{args.directory}/data.zip',ksize=args.ksize,scaled=args.scaled_input,out_csv=f'{args.directory}/results_protein_{args.ksize}.csv',molecule='protein',cores=args.cores, threshold=args.threshold)
            ### Produce csv file with nt and protein containments with FMH OMEGA estimates
            report_dnds = dnds.report_dNdS_pairwise(f"{args.directory}/results_dna_{dna_k}.csv",f"{args.directory}/results_protein_{args.ksize}.csv",ksize=args.ksize)
            report_dnds.to_csv(f'{args.directory}/fmh_omega_{args.ksize}.csv')
        elif args.mode == 'sngl':
            sourmash_ext.sketch_genome_dna(fasta=args.dna_fasta, ksize=dna_k, scaled=args.scaled_input, out_sigfile=f'{args.directory}/dna.zip')
            sourmash_ext.sketch_genome_protein(fasta=args.protein_fasta, ksize=args.ksize, scaled=args.scaled_input, out_sigfile=f'{args.directory}/protein.zip')
            ### Run pairwise instead to estimate cfracs
            sourmash_ext.compare_signatures(ref=f'{args.directory}/dna.zip',query=f'{args.directory}/dna.zip',ksize=dna_k,molecule='dna',working_dir=args.directory)
            sourmash_ext.compare_signatures(ref=f'{args.directory}/protein.zip',query=f'{args.directory}/protein.zip',ksize=args.ksize,molecule='protein',working_dir=args.directory)
            #Report containments
            nt_df = helperfuncs.extract_containment_matrix(mat_csv=f'{args.directory}/compare.dna.{dna_k}.csv')
            #print(nt_df)
            protein_df = helperfuncs.extract_containment_matrix(mat_csv=f'{args.directory}/compare.protein.{args.ksize}.csv')
            #print(protein_df)
            ### Produce csv file with nt and protein containments with FMH OMEGA estimates
            report_df = dnds.report_dNdS(nt_df,protein_df,ksize=args.ksize) #constant is incorrected
            report_df.to_csv(f'{args.directory}/fmh_omega_{args.ksize}.csv') #does not include p_nt_mut nor p_no_mut
            #report_df.to_csv(f'{args.directory}/fmh_omega_{args.ksize}_with_additional_calculations.csv')
        
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
        Use a scale factor of at least 10 for thousands of genomes'
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
        default=0.05,
        nargs='?',
        const='arg_was_not_given',
        help = 'Set containment threshold for sourmash plugin branchwater commands.'
    )    

    parser.add_argument(
        '--mode',
        type=str,
        default="bwmult",
        help="Enables the use of multithreading from sourmash branchwater plugin"
        #help = 'Identify mode to run fmh_omega as sngl, mult, bwmult, bwpair'
    )

''' ############### DEPRECATED ARGUMENTS
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
        '--translate',
        type=str,
        help='indicate yes or no to translate coding sequence'
    )

    parser.add_argument(
        '--outname',
        help = 'Name your study for FMH OMEGA Estimations. Prefix for output files.'
    )
    DEPRECATED ARGUMENTS ###############
'''


    args = parser.parse_args()

    main(args)

