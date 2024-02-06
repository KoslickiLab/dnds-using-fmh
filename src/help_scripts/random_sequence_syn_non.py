import pandas as pd
import argparse
from dNdS.create_random_nt_sequence_using_NG_assumption import get_coding_sequence_from_nucleotide_sequence, positive_selection_based_on_mutation_rate_p, negative_selection_based_on_mutation_rate_p
import random
from Bio.Seq import Seq

WD='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/'
#sequence length
str_len = 10002

#mutation rate p for mutation purposes
mutation_rate_p = 0.01

#output files
RANDOM_REF = open(f'{WD}ref_{str_len}.fna','w')
RANDOM_REF_PROT = open(f'{WD}ref_translated_{str_len}.faa','w')

POSITIVE_SELECTION_QUERIES_NT = open(f'{WD}positive_selection_queries_{str_len}_{mutation_rate_p}.fna','w')
POSITIVE_SELECTION_QUERIES_PROTEIN = open(f'{WD}positive_selection_translated_queries_{str_len}_{mutation_rate_p}.faa','w')

NEGATIVE_SELECTION_QUERIES_NT = open(f'{WD}negative_selection_queries_{str_len}_{mutation_rate_p}.fna','w')
NEGATIVE_SELECTION_QUERIES_PROTEIN = open(f'{WD}negative_selection_translated_queries_{str_len}_{mutation_rate_p}.faa','w')

# Create 100 random mutated sequences
ITERATIONS = 100

#Create a random X nt long sequence for simulation
REF = ''.join(random.choices(['A', 'C', 'G', 'T'], k=str_len))

#Get coding sequence and filter stop codons
REF_coding_sequence = ''.join(get_coding_sequence_from_nucleotide_sequence(REF))
#print('ref: ' + str(len(REF_coding_sequence)))

#Save nt ref sequence to the following files
RANDOM_REF.write(f'>ref_gene\n')
RANDOM_REF.write(f'{REF_coding_sequence}\n')

POSITIVE_SELECTION_QUERIES_NT.write(f'>ref_gene\n')
POSITIVE_SELECTION_QUERIES_NT.write(f'{REF_coding_sequence}\n')

NEGATIVE_SELECTION_QUERIES_NT.write(f'>ref_gene\n')
NEGATIVE_SELECTION_QUERIES_NT.write(f'{REF_coding_sequence}\n')

#save translated ref sequence to the following files
ref_seq_translated = Seq(REF_coding_sequence).translate()
#print('ref_translated: ' + str(len(ref_seq_translated)))
RANDOM_REF_PROT.write(f'>ref_gene\n')
RANDOM_REF_PROT.write(str(ref_seq_translated)+"\n")

for i in range(ITERATIONS):

    #print('ref: ' + str(len(REF_coding_sequence)))
    #print('ref_translated: ' + str(len(ref_seq_translated)))


    #ref sequence is mutated with mutation rate p
    query_positive_nt_seq = positive_selection_based_on_mutation_rate_p(REF_coding_sequence,float(mutation_rate_p))
    POSITIVE_SELECTION_QUERIES_NT.write(f'>positive_{mutation_rate_p}_{i}\n')
    POSITIVE_SELECTION_QUERIES_NT.write(f'{query_positive_nt_seq}\n')

    #Translated mutated sequences 
    translated_positive_queries = Seq(query_positive_nt_seq).translate()
    #print(len(translated_positive_queries))
    POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'>positive_{mutation_rate_p}_{i}\n')
    POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'{translated_positive_queries}\n')

    #ref sequence is mutated with mutation rate p
    query_negative_nt_seq = negative_selection_based_on_mutation_rate_p(REF_coding_sequence,float(mutation_rate_p))
    NEGATIVE_SELECTION_QUERIES_NT.write(f'>negative_{mutation_rate_p}_{i}\n')
    NEGATIVE_SELECTION_QUERIES_NT.write(f'{query_negative_nt_seq}\n')

    #Translated mutated sequences 
    translated_negative_queries = Seq(query_negative_nt_seq).translate()
    #print(len(translated_negative_queries))
    NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'>negative_{mutation_rate_p}_{i}\n')
    NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'{translated_negative_queries}\n')

