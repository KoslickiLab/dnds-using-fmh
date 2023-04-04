import pandas as pd
import argparse
from dNdS import create_ground_truth_file
from Bio.Seq import Seq
from Bio import SeqIO

"""CFrac dN/dS estimator is tested with a ground truth file. This program creates a ground truth file that has\
three columns that include the refernce sequence name, query name, dN, dS, dNdS value between KO genes, and selection"""

def main(args):

    WD = args.wd
    fna_ref_file = args.reference_input
    fna_query_file = args.query_input
    GROUND_TRUTH = args.ground_truth_output

    GROUND_TRUTH_REF_PROT_OUTPUT = open(f'{WD}translated_ref_seq.faa','w')
    GROUND_TRUTH_PROT_QUERIES = open(f'{WD}translated_mutated_queries_seq.faa','w')

    dNdS_report = []

    query_lines=open(fna_query_file,'r').readlines()

    for seq_record in SeqIO.parse(fna_ref_file, "fasta"):
        ref_name = str(seq_record.id)
        ref_seq = str(seq_record.seq)

    for seq_record in SeqIO.parse(fna_query_file, "fasta"):
        query_name = str(seq_record.id)
        query_nt_seq = str(seq_record.seq)

                            
        #get total numver of nucleotide mutations
        total_cdn_mutations = create_ground_truth_file.total_codon_mutations(ref_seq,query_nt_seq) 

        #get total number of nonsynonymous mutations
        total_nonsyn_mutations = create_ground_truth_file.total_aa_differences(ref_seq,query_nt_seq)
        dN = total_nonsyn_mutations/(len(ref_seq)/3)

        #get total number of synonymous mutations
        total_syn_mutations = create_ground_truth_file.total_synonymous_mutations(total_cdn_mutations,total_nonsyn_mutations)
        dS = total_syn_mutations/(len(ref_seq)/3)

        #estimate dNdS using Koslicki's suggestiion
        dNdS = create_ground_truth_file.koslicki_dnds(total_nonsyn_mutations,total_syn_mutations)

        #selection
        if dNdS is None:
            selection='undetermined'
        elif dNdS == 1:
            selection='neutral'
        elif dNdS > 1:
            selection='positive/constrained'
        elif dNdS < 1:
            selection='negative'

        #save report in a list object to later convert into a dataframe
        dNdS_report.append([ref_name, query_name, dN, dS, dNdS, selection])

        #Translated mutated sequences 
        translated_queries = Seq(query_nt_seq).translate()
        GROUND_TRUTH_PROT_QUERIES.write(f'>{query_name}\n')
        GROUND_TRUTH_PROT_QUERIES.write(f'{translated_queries}\n')

        #save translated ref sequence to the following files
        ref_seq_translated = Seq(ref_seq).translate()
        GROUND_TRUTH_REF_PROT_OUTPUT.write(f'>{ref_name}\n')
        GROUND_TRUTH_REF_PROT_OUTPUT.write(str(ref_seq_translated)+"\n")


    #save report to csv output file
    pd.DataFrame(dNdS_report, columns=['ref_gene', 'query_gene','dN','dS','dNdS_ground_truth','selection']).to_csv(f'{WD}{GROUND_TRUTH}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'CFrac dN/dS estimator is tested with a ground truth file. This program creates a ground truth file that has\
        three columns that include the refernce sequence name, mutation rate p, and dNdS value based on the mutation rate p.'
    )

    parser.add_argument(
        '--reference_input',
        type=str,
        help = 'provide fna file of reference sequences to be mutated.'
    )

    parser.add_argument(
        '--query_input',
        type=str,
        help = 'provide fna file of query sequences to compare against the reference'
    )
    
    parser.add_argument(
        '--ground_truth_output',
        type=str,
        help = 'provide the csv filename for ground truth output file.'
    )

    parser.add_argument(
        '--wd',
        type=str,
        default='',
        help = 'provide working directory to save all output'
    )

    args = parser.parse_args()

    main(args)
