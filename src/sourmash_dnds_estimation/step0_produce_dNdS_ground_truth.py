import pandas as pd
import argparse
from dNdS import create_ground_truth_file

"""CFrac dN/dS estimator is tested with a ground truth file. This program creates a ground truth file that has\
three columns that include the refernce sequence name, mutation rate p, and dNdS value based on the mutation rate p"""

def main(args):

    fna_file = args.reference_input
    GROUND_TRUTH = args.ground_truth_output
    mutation_p_list = args.mutation_rate_p.split(',')
    dNdS_report = []
    print(mutation_p_list)

    #loop through nt file to mutate sequences
    with open(fna_file) as file:
        for line in file:
                if line[0] == ">":
                    ref_names = line[1:]
                else:
                    for p in mutation_p_list:
                        #ref sequence is mutated with mutation rate p
                        query_nt_seq = create_ground_truth_file.mutated_sequence_based_on_mutation_rate_p(line,float(p))

                        #get total numver of nucleotide mutations
                        total_nt_mutations = create_ground_truth_file.total_nucleotide_mutations(line,query_nt_seq) 

                        #get total number of nonsynonymous mutations
                        total_nonsyn_mutations = create_ground_truth_file.total_aa_differences(line,query_nt_seq)

                        #get total number of synonymous mutations
                        total_syn_mutations = create_ground_truth_file.total_synonymous_mutations(total_nt_mutations,total_nonsyn_mutations)

                        #estimate dNdS using Koslicki's suggestiion
                        dNdS = create_ground_truth_file.koslicki_dnds(total_nonsyn_mutations,total_syn_mutations,len(line))

                        #save report in a list object to later convert into a dataframe
                        dNdS_report.append([ref_names, p, dNdS])
    file.close()

    #save report to csv output file
    pd.DataFrame(dNdS_report, columns=['ref_sequence_name', 'mutation_rate_p','dNdS_ground_truth']).to_csv(GROUND_TRUTH)

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
        '--mutation_rate_p',
        default='0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1',
        type=str,
        help = 'provide a list of mutation rate p values to mutate reference sequences to obtain ground truth query mutated sequence'
    )

    args = parser.parse_args()

    main(args)
