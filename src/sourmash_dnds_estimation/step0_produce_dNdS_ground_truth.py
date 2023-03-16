import pandas as pd
import argparse
from dNdS import create_ground_truth_file
from Bio.Seq import Seq

"""CFrac dN/dS estimator is tested with a ground truth file. This program creates a ground truth file that has\
three columns that include the refernce sequence name, mutation rate p, and dNdS value based on the mutation rate p"""

def main(args):

    WD = args.wd
    fna_file = args.reference_input
    GROUND_TRUTH = args.ground_truth_output
    GROUND_TRUTH_QUERIES = open(f'{WD}{args.ground_truth_queries_output}','w')
    GROUND_TRUTH_REF_OUTPUT = open(f'{WD}{args.ground_truth_ref_output}','w')
    mutation_p = args.mutation_rate_p
    ITERATIONS = args.iterations

    GROUND_TRUTH_REF_PROT_OUTPUT = open(f'{WD}translated_ref_seq.faa','w')
    GROUND_TRUTH_PROT_QUERIES = open(f'{WD}translated_mutated_queries_seq.faa','w')

    dNdS_report = []

    #loop through nt file to mutate sequences
    #with open(fna_file,'r',encoding='utf-8') as file:
    #GROUND_TRUTH_QUERIES = open('ground_truth_queries.fna','w')
    #GROUND_TRUTH_QUERIES = open('ground_truth_queries.fna','w')

    lines=open(fna_file,'r').readlines()
    for line in lines:
        if line[0] == ">":
            ref_names = line[1:-2] #ignore '>' symbol and '\m'
        else:
            #next funcstion needs ref to be the same size, we ignore the start and stop codons and capitalize all characters        
            ref_seq = line.strip()[3:-3].upper()

            for i in range(ITERATIONS):

                #ref sequence is mutated with mutation rate p
                query_nt_seq = create_ground_truth_file.mutated_sequence_based_on_mutation_rate_p(ref_seq,float(p))

                GROUND_TRUTH_QUERIES.write(f'>gene_{p}_{i}\n')
                GROUND_TRUTH_QUERIES.write(f'{query_nt_seq}\n')

                #Translated mutated sequences 
                translated_queries = Seq(query_nt_seq).translate()
                GROUND_TRUTH_PROT_QUERIES.write(f'>gene_{p}_{i}\n')
                GROUND_TRUTH_PROT_QUERIES.write(f'{translated_queries}\n')
                            
                #get total numver of nucleotide mutations
                total_nt_mutations = create_ground_truth_file.total_nucleotide_mutations(ref_seq,query_nt_seq) 

                #get total number of nonsynonymous mutations
                total_nonsyn_mutations = create_ground_truth_file.total_aa_differences(ref_seq,query_nt_seq)
                dN = total_nonsyn_mutations/(len(ref_seq)/3)

                #get total number of synonymous mutations
                total_syn_mutations = create_ground_truth_file.total_synonymous_mutations(total_nt_mutations,total_nonsyn_mutations)
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
                dNdS_report.append([f'ref_{p}_{i}', p, dN, dS, dNdS, selection, i])

            GROUND_TRUTH_REF_OUTPUT.write(f'>ref_gene\n')
            GROUND_TRUTH_REF_OUTPUT.write(ref_seq+"\n")

            #save translated ref sequence to the following files
            ref_seq_translated = Seq(ref_seq).translate()
            GROUND_TRUTH_REF_PROT_OUTPUT.write(f'>ref_gene\n')
            GROUND_TRUTH_REF_PROT_OUTPUT.write(str(ref_seq_translated)+"\n")


    #save report to csv output file
    pd.DataFrame(dNdS_report, columns=['mutated_sequence', 'mutation_rate_p','dN','dS','dNdS_ground_truth','selection','iteration']).to_csv(GROUND_TRUTH)

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
        '--ground_truth_output',
        type=str,
        help = 'provide the csv filename for ground truth output file.'
    )

    parser.add_argument(
        '--ground_truth_ref_output',
        type=str,
        help = 'generate output files of ref fasta sequences'
    )

    parser.add_argument(
        '--ground_truth_queries_output',
        type=str,
        help = 'generate output files of mutated fasta sequences'
    )

    parser.add_argument(
        '--mutation_rate_p',
        default=0.1,
        type=float,
        help = 'provide a mutation rate p values to mutate reference sequences to obtain ground truth query mutated sequence'
    )

    parser.add_argument(
        '--iterations',
        type=int,
        help = 'provide a list of mutation rate p values to mutate reference sequences to obtain ground truth query mutated sequence'
    )

    parser.add_argument(
        '--wd',
        type=str,
        default='',
        help = 'provide working directory to save all output'
    )

    args = parser.parse_args()

    main(args)
