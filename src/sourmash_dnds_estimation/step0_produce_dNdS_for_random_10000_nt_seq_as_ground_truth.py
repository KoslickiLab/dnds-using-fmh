import pandas as pd
import argparse
from dNdS import create_ground_truth_file
import random
from Bio.Seq import Seq

GROUND_TRUTH = open('/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/dNdS_ground_truth.csv','w')
GROUND_TRUTH_REF_OUTPUT = open('/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/10000_nt_ref_seq.fna','w')
GROUND_TRUTH_QUERIES = open('/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/10000_nt_mutated_queries_seq.fna','w')
GROUND_TRUTH_REF_PROT_OUTPUT = open('/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/10000_prot_ref_seq.fna','w')
GROUND_TRUTH_PROT_QUERIES = open('/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/10000_prot_mutated_queries_seq.fna','w')

ITERATIONS = 100

mutation_p_list = [0.01]

str_len = 10000

REF = ''.join(random.choices(['A', 'C', 'G', 'T'], k=str_len))
ref_seq = REF[3:-3]
GROUND_TRUTH_REF_OUTPUT.write(f'>ref_gene\n')
GROUND_TRUTH_REF_OUTPUT.write(ref_seq+"\n")

ref_seq_translated = Seq(ref_seq).translate()
GROUND_TRUTH_REF_PROT_OUTPUT.write(f'>ref_gene\n')
GROUND_TRUTH_REF_PROT_OUTPUT.write(str(ref_seq_translated)+"\n")

dNdS_report = []

for p in mutation_p_list:

    for i in range(ITERATIONS):

        #ref sequence is mutated with mutation rate p
        query_nt_seq = create_ground_truth_file.mutated_sequence_based_on_mutation_rate_p(REF,float(p))
        GROUND_TRUTH_QUERIES.write(f'>gene_{p}_{i}\n')
        GROUND_TRUTH_QUERIES.write(f'{query_nt_seq}\n')

        translated_queries = Seq(query_nt_seq[3:-3]).translate()
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
        if dNdS == None:
            selection='undetermined'
        elif dNdS == 1:
            selection='neutral'
        elif dNdS > 1:
            selection='positive/constrained'
        elif dNdS < 1:
            selection='negative'

        #save report in a list object to later convert into a dataframe
        dNdS_report.append(['ref_random_10000_nt', p, i, dN, dS, dNdS, selection])

#save report to csv output file
pd.DataFrame(dNdS_report, columns=['ref_sequence_name', 'mutation_rate_p','iteration','dN','dS','dNdS_ground_truth','selection']).to_csv(GROUND_TRUTH)

