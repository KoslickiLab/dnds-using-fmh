#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse
from fmh_omega import helperfuncs,dnds,sourmash_api
import subprocess

def main(args):

    ###Arguments
    dna_fasta = args.cds_input
    klst = args.klist
    s = args.scaled_input
    on = args.outname
    wd = args.working_dir
    single =args.singleton

    #Create signature directory
    subprocess.run(f'mkdir {wd}/signatures', shell=True, check=True)

    #kmers list for sourmash
    sm_dna_klst=helperfuncs.return_dna_klist_parameters(kmer_list=klst)
    sm_protein_klst=helperfuncs.return_protein_klist_parameters(kmer_list=klst)

    #Translate before sketching protein
    helperfuncs.translate_CDS(cds_fasta=dna_fasta, out_name=on, working_dir=wd)
    
    #Sketch signatures
    if single == 'no':
        sourmash_api.sketch_genome_dna(fasta=dna_fasta, klist=sm_dna_klst, scaled=s, out_sigfile=f'{wd}/signatures/{on}.dna.sig.gzip')
        sourmash_api.sketch_genome_protein(fasta=f'{wd}/{on}.translated.fasta', klist=sm_protein_klst, scaled=s, out_sigfile=f'{wd}/signatures/{on}.protein.sig.gzip')
    elif single == 'yes':    
        sourmash_api.sketch_singleton_genome_dna(fasta=dna_fasta, klist=sm_dna_klst, scaled=s, out_sigfile=f'{wd}/signatures/{on}.dna.sig.gzip')
        sourmash_api.sketch_singleton_genome_protein(fasta=f'{wd}/{on}.translated.fasta', klist=sm_protein_klst, scaled=s, out_sigfile=f'{wd}/signatures/{on}.protein.sig.gzip')
    #Concatenate files for single signature
    signature_list=helperfuncs.return_signature_list(working_dir=f'{wd}', molecule='dna')
    if len(signature_list) > 1:
        sourmash_api.cat_signatures(working_dir=f'{wd}', out_name=on, molecule='dna')
        sourmash_api.cat_signatures(working_dir=f'{wd}', out_name=on, molecule='protein')

        #Create compare directory
        subprocess.run(f'mkdir {wd}/compare_dna', shell=True, check=True)
        subprocess.run(f'mkdir {wd}/compare_protein', shell=True, check=True)
        #Containments using sourmash compare
        kmer_list=klst.split(',')
        print(kmer_list)
        for k in kmer_list:
            dna_k = int(k)*3
            sourmash_api.compare_signatures(ref=f'{wd}/signatures/{on}.cat.dna.sig.gzip', query=f'{wd}/signatures/{on}.cat.dna.sig.gzip', ksize=dna_k, molecule='dna', working_dir=f'{wd}')
            sourmash_api.compare_signatures(ref=f'{wd}/signatures/{on}.cat.prot.sig.gzip', query=f'{wd}/signatures/{on}.cat.prot.sig.gzip', ksize=k, molecule='protein', working_dir=f'{wd}')
            #Report containments
            nt_df = helperfuncs.containments(mat_df=f'{wd}/compare_dna/compare.dna.{dna_k}.csv',ksize=k)
            nt_df.to_csv(f'{wd}/nt_containment{dna_k}.csv')
            protein_df = helperfuncs.containments(mat_df=f'{wd}/compare_dna/compare.protein.{k}.csv',ksize=k)
            protein_df.to_csv(f'{wd}/prot_containment{k}.csv')
            ### Produce csv file with nt and protein containments with FMH OMEGA estimates
            report_df = dnds.report_dNdS(nt_df,protein_df)
            report_df.to_csv(f'{wd}/fmh_omega_{k}.csv')
    elif len(signature_list) == 1:
        #Create compare directory
        subprocess.run(f'mkdir {wd}/compare_dna', shell=True, check=True)
        subprocess.run(f'mkdir {wd}/compare_protein', shell=True, check=True)
        #Containments using sourmash compare
        kmer_list=klst.split(',')
        print(kmer_list)
        for k in kmer_list:
            dna_k = int(k)*3
            sourmash_api.compare_signatures(ref=f'{wd}/signatures/{on}.dna.sig.gzip', query=f'{wd}/signatures/{on}.dna.sig.gzip', ksize=dna_k, molecule='dna', working_dir=f'{wd}')
            sourmash_api.compare_signatures(ref=f'{wd}/signatures/{on}.protein.sig.gzip', query=f'{wd}/signatures/{on}.protein.sig.gzip', ksize=k, molecule='protein', working_dir=f'{wd}')
            #Report containments
            nt_df = helperfuncs.containments(mat_df=f'{wd}/compare_dna/compare.dna.{dna_k}.csv',ksize=int(k))
            nt_df.to_csv(f'{wd}/nt_containment{dna_k}.csv')
            protein_df = helperfuncs.containments(mat_df=f'{wd}/compare_protein/compare.protein.{k}.csv',ksize=int(k))
            protein_df.to_csv(f'{wd}/prot_containment{k}.csv')
            ### Produce csv file with nt and protein containments with FMH OMEGA estimates
            report_df = dnds.report_dNdS(nt_df,protein_df)
            report_df.to_csv(f'{wd}/fmh_omega_{k}.csv')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'dN/dS estimator for metagenomic data using the containment index between k-mer sets of genomic samples'
    )

    parser.add_argument(
        '--cds_input',
        help = 'Input of CDS fasta file.\
        This file is a pairwise matrix produced from sourmash compare that includes containment indexes between nucleotide sequences.'
    )

    parser.add_argument(
        '--klist',
        type=str,
        help = 'Identify a list of ksizes used to produce sketches.\
        ksize is required to be the same used in both the containment indexes calculated for nucleotide and protein sequences.'
    )

    parser.add_argument(
        '--scaled_input',
        type=int,
        help = 'Identify a scaled factor for signature sketches.\
        E.g., A scaled factor = 1 will include all k-mers in final signature.'
    )

    parser.add_argument(
        '--singleton',
        type=str,
        help = 'Identify whether you will sketch with singleton option. Yes or No.\
        E.g., yes, when a single fasta file with multiple sequences will be used for analysis.'
    )

    parser.add_argument(
        '--outname',
        help = 'Name your study for FMH OMEGA Estimations. Prefix for output files.'
    )

    parser.add_argument(
        '--working_dir',
        help = 'Output directory for FMH Omega estimation.'
    )

    args = parser.parse_args()

    main(args)
