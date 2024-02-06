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

    #Create signature directory
    subprocess.run(f'mkdir {wd}/signatures', shell=True, check=True)
    #Sketch signatures
    sourmash_api.sketch_dna(fasta=dna_fasta, klist=klst, scaled=s, out_sigfile=f'{wd}/signatures/{on}.dna.sig.gzip')
    #Translate before sketching protein
    helperfuncs.translate_CDS(cds_fasta=dna_fasta, out_name=on)
    sourmash_api.sketch_protein(fasta=f'{on}.translated.fasta', klist=klst, scaled=s, out_sigfile=f'{wd}/signatures/{on}.prot.sig.gzip')
    #Concatenate files for single signature
    sourmash_api.cat_signatures(working_dir=wd, out_name=on, molecule='dna')
    sourmash_api.cat_signatures(working_dir=wd, out_name=on, molecule='protein')

    #Create compare directory
    subprocess.run(f'mkdir {wd}/compare_dna', shell=True, check=True)
    subprocess.run(f'mkdir {wd}/compare_protein', shell=True, check=True)
    #Containments using sourmash compare
    for k in klst:
        sourmash_api.compare_signatures(ref=f'{wd}/signatures/{on}.cat.dna.sig.gzip', query=f'{wd}/signatures/{on}.dna.sig.gzip', ksize=k, molecule='dna', working_dir=wd)
        sourmash_api.compare_signatures(ref=f'{wd}/signatures/{on}.cat.prot.sig.gzip', query=f'{wd}/signatures/{on}.cat.prot.sig.gzip', ksize=k, molecule='protein', working_dir=wd)
        #Report containments
        nt_df = helperfuncs.containments(mat_df=f'{wd}/compare_dna/compare.dna.${k}.csv',ksize=k)
        nt_df.to_csv(f'{wd}/nt_containment{k}.csv')
        protein_df = helperfuncs.containments(mat_df=f'{wd}/compare_dna/compare.protein.${k}.csv',ksize=k)
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
        type=list,
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
        '--outname',
        help = 'Name your study for FMH OMEGA Estimations. Prefix for output files.'
    )

    parser.add_argument(
        '--wd',
        help = 'Output directory for FMH Omega estimation.'
    )

    args = parser.parse_args()

    main(args)
