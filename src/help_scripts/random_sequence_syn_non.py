import pandas as pd
import argparse
import dNdS.create_random_nt_sequence_using_NG_assumption 
import random
from Bio.Seq import Seq

WD='/data/jzr5814/sourmash_dnds_estimation/tests/test/'
#sequence length
str_len = 10002

#mutation rate p for mutation purposes
mutation_p_list = 0.001

#output files
RANDOM_REF = open(f'{WD}ref_{str_len}.fna','w')
RANDOM_REF_PROT = open(f'{WD}ref_translated_{str_len}.faa','w')

POSITIVE_SELECTION_QUERIES_NT = open(f'{WD}positive_selection_queries_{str_len}_{mutation_p_list}.fna','w')
POSITIVE_SELECTION_QUERIES_PROTEIN = open(f'{WD}positive_selection_translated_queries_{str_len}_{mutation_p_list}.faa','w')

NEGATIVE_SELECTION_QUERIES_NT = open(f'{WD}negative_selection_queries_{str_len}_{mutation_p_list}.fna','w')
NEGATIVE_SELECTION_QUERIES_PROTEIN = open(f'{WD}negative_selection_translated_queries_{str_len}_{mutation_p_list}.faa','w')

# Create 100 random mutated sequences
ITERATIONS = 100

#Create a random X nt long sequence for simulation
REF = ''.join(random.choices(['A', 'C', 'G', 'T'], k=str_len))

#Get coding sequence and filter stop codons
REF_coding_sequence = get_coding_sequence_from_nucleotide_sequence(REF)

#Save nt ref sequence to the following files
RANDOM_REF.write(f'>ref_gene\n')
RANDOM_REF.write(REF_coding_sequence+"\n")

#save translated ref sequence to the following files
ref_seq_translated = Seq(REF_coding_sequence).translate()
RANDOM_REF_PROT.write(f'>ref_gene\n')
RANDOM_REF_PROT.write(str(ref_seq_translated)+"\n")

#an empty list to add dN/dS values
dNdS_report = []

for i in range(ITERATIONS):

    #ref sequence is mutated with mutation rate p
    query_nt_seq = positive_selection_based_on_mutation_rate_p(ref_seq,float(p))
    POSITIVE_SELECTION_QUERIES_NT.write(f'>gene_{p}_{i}\n')
    POSITIVE_SELECTION_QUERIES_NT.write(f'{query_nt_seq}\n')

    #Translated mutated sequences 
    translated_queries = Seq(query_nt_seq).translate()
    POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'>gene_{p}_{i}\n')
    POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'{translated_queries}\n')

    #ref sequence is mutated with mutation rate p
    query_nt_seq = positive_selection_based_on_mutation_rate_p(ref_seq,float(p))
    NEGATIVE_SELECTION_QUERIES_NT.write(f'>gene_{p}_{i}\n')
    NEGATIVE_SELECTION_QUERIES_NT.write(f'{query_nt_seq}\n')

    #Translated mutated sequences 
    translated_queries = Seq(query_nt_seq).translate()
    NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'>gene_{p}_{i}\n')
    NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'{translated_queries}\n')

