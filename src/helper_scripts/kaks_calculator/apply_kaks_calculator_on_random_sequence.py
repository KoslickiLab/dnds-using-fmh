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

def main(args):

    input_file = args.input
    input_filename = os.path.basename(args.input)
    wd=args.wd.rstrip('/')
    axt_file = open(f'{wd}/{input_filename}.axt','w')
    
    seguid_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"), lambda rec: seguid(rec.seq))

    for key in seguid_dict:
        ref = seguid_dict[key].description
        ref_seq = str(seguid_dict[key].seq)
        for key in seguid_dict:
            key_name = str(seguid_dict[key].description)
            print(f'{ref}_vs_{key_name}')
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

    args = parser.parse_args()

    main(args)