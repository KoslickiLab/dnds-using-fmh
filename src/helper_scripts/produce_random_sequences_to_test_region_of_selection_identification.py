import pandas as pd
import argparse
from fmh_omega.helperfuncs import get_coding_sequence_from_nucleotide_sequence, positive_selection_based_on_mutation_rate_p, negative_selection_based_on_mutation_rate_p
import random
from Bio.Seq import Seq

WD='/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_practical_considerations/lengths/'
#sequence length
str_len = 10002
shared_len = 1000

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
ITERATIONS = 3

#Create a random X nt long sequence for simulation
REF = ''.join(random.choices(['A', 'C', 'G', 'T'], k=str_len))

#Get coding sequence and filter stop codons
REF_coding_sequence = ''.join(get_coding_sequence_from_nucleotide_sequence(REF))

#Save nt ref sequence to the following files
RANDOM_REF.write(f'>ref_gene\n')
RANDOM_REF.write(f'{REF_coding_sequence}\n')

POSITIVE_SELECTION_QUERIES_NT.write(f'>ref_gene\n')
POSITIVE_SELECTION_QUERIES_NT.write(f'{REF_coding_sequence}\n')

NEGATIVE_SELECTION_QUERIES_NT.write(f'>ref_gene\n')
NEGATIVE_SELECTION_QUERIES_NT.write(f'{REF_coding_sequence}\n')

#save translated ref sequence to the following files
ref_seq_translated = Seq(REF_coding_sequence).translate()
RANDOM_REF_PROT.write(f'>ref_gene\n')
RANDOM_REF_PROT.write(str(ref_seq_translated)+"\n")

POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'>ref_gene\n')
POSITIVE_SELECTION_QUERIES_PROTEIN.write(str(ref_seq_translated)+"\n")

NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'>ref_gene\n')
NEGATIVE_SELECTION_QUERIES_PROTEIN.write(str(ref_seq_translated)+"\n")

for i in range(ITERATIONS):
    #ref sequence is mutated with mutation rate p
    query_positive_nt_seq = positive_selection_based_on_mutation_rate_p(REF_coding_sequence,float(mutation_rate_p))
    POSITIVE_SELECTION_QUERIES_NT.write(f'>positive_{mutation_rate_p}_{i}\n')
    POSITIVE_SELECTION_QUERIES_NT.write(f'{query_positive_nt_seq}\n')

    #Translated mutated sequenceså
    translated_positive_queries = Seq(query_positive_nt_seq).translate()
    POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'>positive_{mutation_rate_p}_{i}\n')
    POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'{translated_positive_queries}\n')

    #ref sequence is mutated with mutation rate p
    query_negative_nt_seq = negative_selection_based_on_mutation_rate_p(REF_coding_sequence,float(mutation_rate_p))
    NEGATIVE_SELECTION_QUERIES_NT.write(f'>negative_{mutation_rate_p}_{i}\n')
    NEGATIVE_SELECTION_QUERIES_NT.write(f'{query_negative_nt_seq}\n')

    #Translated mutated sequences 
    translated_negative_queries = Seq(query_negative_nt_seq).translate()
    NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'>negative_{mutation_rate_p}_{i}\n')
    NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'{translated_negative_queries}\n')

    #Create a random region X nt long that is shared in both selected sequences
    shared_len = random.randint(1000, 5000)
    SHARED_REGION = ''.join(random.choices(['A','C','G','T'], k=shared_len))

    #insert shared region in positive selected sequence
    positive_position = random.randint(2500, 7500)
    POSITIVE_SHARED_SEQUENCE = query_positive_nt_seq.replace(SHARED_REGION[positive_position:len(SHARED_REGION)],SHARED_REGION)
    POSITIVE_SELECTION_QUERIES_NT.write(f'>positive_shared_{mutation_rate_p}_{i}\n')
    POSITIVE_SELECTION_QUERIES_NT.write(f'{query_positive_nt_seq}\n')

    #Translated with shared region 
    translated_positive_shared_queries = Seq(POSITIVE_SHARED_SEQUENCE).translate()
    POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'>positive_shared_{mutation_rate_p}_{i}\n')
    POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'{translated_positive_shared_queries}\n')

    #Insert shared region in negative selected sequence
    negative_position = random.randint(2500, 7500)
    NEGATIVE_SHARED_SEQUENCE = query_negative_nt_seq.replace(SHARED_REGION[negative_position:len(SHARED_REGION)],SHARED_REGION)
    NEGATIVE_SELECTION_QUERIES_NT.write(f'>negative_shared_{mutation_rate_p}_{i}\n')
    NEGATIVE_SELECTION_QUERIES_NT.write(f'{query_negative_nt_seq}\n')

    #Translated with shared region 
    translated_negative_shared_queries = Seq(NEGATIVE_SHARED_SEQUENCE).translate()
    NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'>negative_shared_{mutation_rate_p}_{i}\n')
    NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'{translated_negative_shared_queries}\n')




