import pandas as pd
import argparse
from helper_functions.helperfuncs import get_coding_sequence_from_nucleotide_sequence, positive_selection_based_on_mutation_rate_p, negative_selection_based_on_mutation_rate_p
import random
from Bio.Seq import Seq


def main(args):
    #output files
    RANDOM_REF = open(f'{args.wd.rstrip('/')}/ref_{args.length}.fna','w')
    RANDOM_REF_PROT = open(f'{args.wd.rstrip('/')}/ref_translated_{args.length}.faa','w')

    POSITIVE_SELECTION_QUERIES_NT = open(f'{args.wd.rstrip('/')}/positive_selection_queries_{args.length}_{args.prate}.fna','w')
    POSITIVE_SELECTION_QUERIES_PROTEIN = open(f'{args.wd.rstrip('/')}/positive_selection_translated_queries_{args.length}_{args.prate}.faa','w')

    NEGATIVE_SELECTION_QUERIES_NT = open(f'{args.wd.rstrip('/')}/negative_selection_queries_{args.length}_{args.prate}.fna','w')
    NEGATIVE_SELECTION_QUERIES_PROTEIN = open(f'{args.wd.rstrip('/')}/negative_selection_translated_queries_{args.length}_{args.prate}.faa','w')

    #Create a random X nt long sequence for simulation
    REF = ''.join(random.choices(['A', 'C', 'G', 'T'], k=args.length))

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

    for i in range(args.iterations):
        #ref sequence is mutated with mutation rate p
        query_positive_nt_seq = positive_selection_based_on_mutation_rate_p(REF_coding_sequence,float(args.prate))
        POSITIVE_SELECTION_QUERIES_NT.write(f'>positive_{args.prate}_{i}\n')
        POSITIVE_SELECTION_QUERIES_NT.write(f'{query_positive_nt_seq}\n')

        #Translated mutated sequences 
        translated_positive_queries = Seq(query_positive_nt_seq).translate()
        POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'>positive_{args.prate}_{i}\n')
        POSITIVE_SELECTION_QUERIES_PROTEIN.write(f'{translated_positive_queries}\n')

        #ref sequence is mutated with mutation rate p
        query_negative_nt_seq = negative_selection_based_on_mutation_rate_p(REF_coding_sequence,float(args.prate))
        NEGATIVE_SELECTION_QUERIES_NT.write(f'>negative_{args.prate}_{i}\n')
        NEGATIVE_SELECTION_QUERIES_NT.write(f'{query_negative_nt_seq}\n')

        #Translated mutated sequences 
        translated_negative_queries = Seq(query_negative_nt_seq).translate()
        NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'>negative_{args.prate}_{i}\n')
        NEGATIVE_SELECTION_QUERIES_PROTEIN.write(f'{translated_negative_queries}\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Script to generate selection simulation for random sequences of a specific length and prate'
    )

    parser.add_argument(
        '--wd',
        nargs='?',
        const='arg_was_not_given',
        type=str,
        help = 'Indicate working directory to save newly generated random sequences'
    )

    parser.add_argument(
        '--length',
        nargs='?',
        const='arg_was_not_given',
        type=int,
        help = 'length of interest'
    )

    parser.add_argument(
        '--prate',
        nargs='?',
        const='arg_was_not_given',
        type=str,
        help = 'mutation rate for simulation'
    )

    parser.add_argument(
        '--iterations',
        nargs='?',
        const='arg_was_not_given',
        type=int,
        help = 'total sequences to be generated at prate (has been ran 100 times)'
    )

    args = parser.parse_args()

    main(args)


