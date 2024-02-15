#!/usr/bin/env python3
"""New approach to estimaating dN/dS ratio of metagenomic data"""

import argparse
from fmh_omega import helperfuncs,dnds,sourmash_ext
import subprocess

def main(args):

    ###Arguments
    dna_fasta = args.cds_input_list
    klst = args.klist
    s = args.scaled_input
    on = args.outname
    wd = args.working_dir
    multiple =args.multiple

    #Create signature directory
    subprocess.run(f'mkdir {wd}/signatures', shell=True, check=True)
    #kmers list for sourmash
    sm_dna_klst=helperfuncs.return_dna_klist_parameters(kmer_list=klst)
    sm_protein_klst=helperfuncs.return_protein_klist_parameters(kmer_list=klst)
    
    #store file information in lists
    fastn_files=[] #dna fasta file
    fasta_files=[] #protein fasta files
    #fast_names = [] #names of entries
    for files in [line.strip() for line in open(f'{dna_fasta}', 'r')]:
        fastn_files.append(files.split(',')[0])
        #fast_names.append(files.split(',')[1])
        fasta_files.append(f'{wd}/{files.split(',')[0].split(".")[0]}.translated.fasta')
        
    #Translate before sketching protein
    for pos in range(len(fastn_files)):
        helperfuncs.translate_CDS(cds_fasta=f'{wd}/{fastn_files[pos]}', out_name=f'{fasta_files[pos]}')

    #Sketch signatures
    fastn_lst = f'{wd}/'+f' {wd}/'.join(fastn_files)
    fasta_lst = f' '.join(fasta_files)
    dna_sig_outname = f'{wd}/signatures/{on}.dna.sig.gzip'
    prot_sig_outname = f'{wd}/signatures/{on}.protein.sig.gzip'
    if multiple == 'yes': #input are multiple fastn files
        sourmash.sketch_genome_dna(fasta=fastn_lst, klist=sm_dna_klst, scaled=s, out_sigfile=dna_sig_outname, multiple={multiple})
        sourmash.sketch_genome_protein(fasta=fasta_lst, klist=sm_protein_klst, scaled=s, out_sigfile=prot_sig_outname, multiple={multiple})
    elif multiple == 'no':
        sourmash.sketch_genome_dna(fasta=fastn_lst, klist=sm_dna_klst, scaled=s, out_sigfile=dna_sig_outname)
        sourmash.sketch_genome_protein(fasta=fasta_lst, klist=sm_protein_klst, scaled=s, out_sigfile=prot_sig_outname)
    
    #Create compare directory
    subprocess.run(f'mkdir {wd}/compare_dna', shell=True, check=True)
    subprocess.run(f'mkdir {wd}/compare_protein', shell=True, check=True)
    
    #Containments using sourmash compare
    kmer_list=klst.split(',')
    for k in kmer_list:
        dna_k = int(k)*3
        sourmash.compare_signatures(ref=f'{wd}/signatures/{on}.dna.sig.gzip', query=f'{wd}/signatures/{on}.dna.sig.gzip', ksize=dna_k, molecule='dna', working_dir=f'{wd}')
        sourmash.compare_signatures(ref=f'{wd}/signatures/{on}.protein.sig.gzip', query=f'{wd}/signatures/{on}.protein.sig.gzip', ksize=k, molecule='protein', working_dir=f'{wd}')
        #Report containments
        nt_df = helperfuncs.containments(mat_df=f'{wd}/compare_dna/compare.dna.{dna_k}.csv',ksize=int(k),multiple=multiple)
        protein_df = helperfuncs.containments(mat_df=f'{wd}/compare_protein/compare.protein.{k}.csv',ksize=int(k),multiple=multiple)
        ### Produce csv file with nt and protein containments with FMH OMEGA estimates
        report_df = dnds.report_dNdS(nt_df,protein_df)
        report_df.to_csv(f'{wd}/fmh_omega_{k}.csv')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'dN/dS estimator for metagenomic data using the containment index between k-mer sets of genomic samples'
    )

    parser.add_argument(
        '--cds_input_list',
        help = 'Input txt file that contains fasta files for sketching.\
        In first column, have fasta file name and second column have nickname for each fasta files.'
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
        '--multiple',
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

