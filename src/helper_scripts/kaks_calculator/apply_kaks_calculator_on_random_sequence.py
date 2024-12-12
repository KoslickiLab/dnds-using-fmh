#!/usr/bin/env python3

'''
    This script is to create a AXT file from a random sample of sequences to run pairwise dn/ds estimates using the kaks_calculator. 
    If you have a fasta file with multiple entries, then you want to use this script to create an axt file for pairwise dn/ds.
    This script was mean't to run on a fasta file generated with a random sequence of length X. The length has to be divisible by 3 or the kaks_calculator will return an error.
    If you have sequences of differing lengths, then using an aligner such as clustalw can help but another script would have to be implemented to create a AXT file for kaks_calculator input. 
'''

import argparse
import subprocess
import os
from Bio.SeqUtils.CheckSum import seguid
from Bio import SeqIO
from collections import defaultdict

def main(args):

    input_file = args.input
    input_filename = os.path.basename(args.input)
    wd=args.wd.rstrip('/')
    axt_file = open(f'{wd}/{input_filename}.axt','w')

    seq_dict = defaultdict(list)
    for record in SeqIO.parse(input_file, "fasta"):
        key = seguid(record.seq)
        seq_dict[key].append(record)

    if args.short:
        # Shorter sequences create duplicate keys when running SeqIO.to_dict, Assuming seq_dict contains lists of SeqRecord objects
        for key, records in seq_dict.items():  # Iterate through each key and list of SeqRecord objects
            for ref_record in records:  # Loop over each SeqRecord in the list of the current key
                ref = ref_record.description  # Access the description of the current reference record
                ref_seq = str(ref_record.seq)  # Access the sequence of the current reference record
                
                # Now compare this reference record with every other SeqRecord in seq_dict
                for other_key, other_records in seq_dict.items():  # Compare with every key in the dictionary
                    for other_record in other_records:  # Loop through the list of SeqRecord objects for the other key
                        key_name = str(other_record.description)  # Get the description of the other SeqRecord
                        
                        # Print and write comparison information
                        axt_file.write(f'{ref}_vs_{key_name}\n')  # Write the comparison to the file
                        axt_file.write(ref_seq + '\n')  # Write the reference sequence to the file
                        axt_file.write(str(other_record.seq) + '\n')  # Write the other sequence to the file
                        axt_file.write('\n')  # Write an empty line to separate comparison

    else:
        seguid_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"), lambda rec: seguid(rec.seq))
        for key in seguid_dict:
            ref = seguid_dict[key].description
            ref_seq = str(seguid_dict[key].seq)
            for key in seguid_dict:
                key_name = str(seguid_dict[key].description)
                axt_file.write(f'{ref}_vs_{key_name}')
                axt_file.write('\n')
                axt_file.write(ref_seq)
                axt_file.write('\n')
                axt_file.write(str(seguid_dict[key].seq))
                axt_file.write('\n')
                axt_file.write('\n')

    cmd=f'KaKs_Calculator -i {wd}/{input_filename}.axt -o {wd}/{input_filename}.axt.kaks -m {args.method}'
    subprocess.run(cmd, shell=True, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Run KaKs_calculator on a random fasta file generated and simulated for testing ground truths'
    )

    parser.add_argument(
        '--input',
        type=str,
        help = 'Input fasta file of randomly generated sequences, make sure your reference is the first entry.'
    )

    parser.add_argument(
        '--method',
        type=str,
        help = 'Choose one from the following methods: GNG, GY, LPB, LWL, NG, YN.'
    )

    parser.add_argument(
        '--wd',
        type=str,
        help = 'Identify the working directory where output files will be produced.'
    )
    
    parser.add_argument(
        '--short',
        action='store_true',  # This makes the argument a boolean flag
        help='When applying on shorter sequences, less than 1000 nucleotides long (default: False)'
    )

    args = parser.parse_args()

    main(args)