#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse
from fmh_omega import helperfuncs,dnds,sourmash_ext
import subprocess
import time
import multiprocessing as mp
import numpy as np
import math

def main(args):
    
    ### ARGUMENTS
    dna_fasta = args.fasta_input_list
    k = args.ksize
    s = args.scaled_input
    on = args.outname
    wd = args.directory
    m =args.mode
    translate_cds=args.translate
    total_cores =args.cores

    ### PREPARING RUN
    #kmers list for sourmash
    #sm_dna_klst=helperfuncs.return_dna_klist_parameters(kmer_list=klst)
    #sm_protein_klst=helperfuncs.return_protein_klist_parameters(kmer_list=klst)
    #Containments using sourmash compare and multisearch
    #kmer_list=klst.split(',')
    dna_k = k*3
    #get total expected signatures
    total_num_signatures=-1
    with open(f'{dna_fasta}') as infp:
        for line in infp:
            if line.strip():
                total_num_signatures += 1
    #total cores to use
    #total_num_signatures = total_num_signatures-1 #subtract 1 for the header
    #total_cores = min(mp.cpu_count(), math.ceil(total_num_signatures/1000))
    #total_cores = 100 #testing

    ### RUN WHEN NOT USING SOURMASH BRANCHWATER PLUGIN
    if m != "bwmult" and m !="bwpair":
        #store file information in lists
        fastn_files=[] #dna fasta file
        fasta_files=[] #protein fasta files
        for files in [line.strip() for line in open(f'{dna_fasta}', 'r')]:
            if files != '' and 'genome_filename' not in files:
                fastn_filename = files.split(',')[1]
                name = files.split(',')[1].split('.')[0]
                fastn_files.append(fastn_filename)
                fasta_files.append(f'{name}.translated.fasta')
        if translate_cds == 'yes':
        #Translate before sketching protein
            for pos in range(len(fastn_files)):
                helperfuncs.translate_CDS(cds_fasta=f'{fastn_files[pos]}', out_name=f'{fasta_files[pos]}')

        #Create signature directory
        subprocess.run(f'mkdir {wd}/signatures', shell=True, check=True)
        #Sketch signatures
        fastn_lst = f'{wd}/'+f' {wd}/'.join(fastn_files)
        fasta_lst = f' '.join(fasta_files)
        dna_sig_outname = f'{wd}/signatures/{on}.dna.sig.gzip'
        prot_sig_outname = f'{wd}/signatures/{on}.protein.sig.gzip'
        if m == 'mult': #input are fastn files
            sourmash_ext.sketch_genome_dna(fasta=fastn_lst, ksize=dna_k, scaled=s, out_sigfile=dna_sig_outname, multiple={m})
            sourmash_ext.sketch_genome_protein(fasta=fasta_lst, ksize=k, scaled=s, out_sigfile=prot_sig_outname, multiple={m})
        elif m == 'sngl':
            sourmash_ext.sketch_genome_dna(fasta=fastn_lst, ksize=dna_k, scaled=s, out_sigfile=dna_sig_outname)
            sourmash_ext.sketch_genome_protein(fasta=fasta_lst, ksize=k, scaled=s, out_sigfile=prot_sig_outname)
        # Create compare directory
        #subprocess.run(f'mkdir {wd}/compare_dna', shell=True, check=True)
        #subprocess.run(f'mkdir {wd}/compare_protein', shell=True, check=True)
        ### run sourmash comapre for cfracs
        sourmash_ext.compare_signatures(ref=f'{wd}/signatures/{on}.dna.sig.gzip', query=f'{wd}/signatures/{on}.dna.sig.gzip', ksize=dna_k, molecule='dna', working_dir=f'{wd}')
        sourmash_ext.compare_signatures(ref=f'{wd}/signatures/{on}.protein.sig.gzip', query=f'{wd}/signatures/{on}.protein.sig.gzip', ksize=k, molecule='protein', working_dir=f'{wd}')
        #Report containments
        nt_df = helperfuncs.containments(mat_df=f'{wd}/compare.dna.{dna_k}.csv',ksize=int(k),multiple=m)
        protein_df = helperfuncs.containments(mat_df=f'{wd}/compare.protein.{k}.csv',ksize=int(k),multiple=m)
        ### Produce csv file with nt and protein containments with FMH OMEGA estimates
        report_df = dnds.report_dNdS(nt_df,protein_df)
        report_df.to_csv(f'{wd}/fmh_omega_{k}.csv')

    ### RUN WHEN USING SOURMASH BRANCHWATER PLUGIN
    elif m == "bwmult" or m == "bwpair":
        #sourmash_ext.run_manysketch(fasta_file_csv=dna_fasta, ksize=dna_k, scaled=s, molecule='dna',cores=total_cores, working_dir=wd)
        #sourmash_ext.run_manysketch(fasta_file_csv=dna_fasta, ksize=k, scaled=s, molecule='protein',cores=total_cores, working_dir=wd)
        sourmash_ext.run_manysketch(fasta_file_csv=dna_fasta, ksize=k, scaled=s, cores=total_cores, working_dir=wd)
        if m == "bwmult":
            ### Run multisearch to estimate cfracs
            sourmash_ext.run_multisearch(ref_zipfile=f'{wd}/dna.zip',query_zipfile=f'{wd}/dna.zip',ksize=dna_k,scaled=s,out_csv=f'{wd}/results_dna_{dna_k}.csv',molecule='DNA',cores=total_cores)
            sourmash_ext.run_multisearch(ref_zipfile=f'{wd}/protein.zip',query_zipfile=f'{wd}/protein.zip',ksize=k,scaled=s,out_csv=f'{wd}/results_protein_{k}.csv',molecule='protein',cores=total_cores)
            ### Produce csv file with nt and protein containments with FMH OMEGA estimates
            report_dnds = dnds.report_dNdS_multisearch(f"{wd}/results_dna_{dna_k}.csv",f"{wd}/results_protein_{k}.csv",ksize=k)
            report_dnds.to_csv(f'{wd}/fmh_omega_{k}.csv')
        elif m == "bwpair":
            ### Run pairwise instead to estimate cfracs
            #sourmash_ext.run_pairwise(zipfile=f'{wd}/dna.zip',ksize=dna_k,scaled=s,out_csv=f'{wd}/results_dna_{dna_k}.csv',molecule='DNA',cores=total_cores)
            #sourmash_ext.run_pairwise(zipfile=f'{wd}/protein.zip',ksize=k,scaled=s,out_csv=f'{wd}/results_protein_{k}.csv',molecule='protein',cores=total_cores)
            sourmash_ext.run_pairwise(zipfile=f'{wd}/data.zip',ksize=k,scaled=s,out_csv=f'{wd}/results_protein_{k}.csv',molecule='protein',cores=total_cores)
            ### Produce csv file with nt and protein containments with FMH OMEGA estimates
            report_dnds = dnds.report_dNdS_pairwise(f"{wd}/results_dna_{dna_k}.csv",f"{wd}/results_protein_{k}.csv",ksize=k)
            report_dnds.to_csv(f'{wd}/fmh_omega_{k}.csv')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'dN/dS estimator for metagenomic data using the containment index between k-mer sets of genomic samples'
    )

    parser.add_argument(
        '--fasta_input_list',
        help = 'Input txt file that contains fasta files for sketching.\
        In first column, have fasta file name and second column have nickname for each fasta files.'
    )

    parser.add_argument(
        '--ksize',
        type=int,
        help = 'Identify a ksize used to produce sketches.\
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

    args = parser.parse_args()

    main(args)

