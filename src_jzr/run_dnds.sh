#!/bin/bash
set -eoux pipefail

#Test run on small dataset
#nohup python3 dnds_using_containment_index.py --ref_fasta ../../data/tests/50_sequences.fna --query_fasta ../../data/tests/5_nt_sequences.fna --predict frame --k 7,14,21,28 --wd ../../data/tests/ > log.txt 2>&1 &

#Running sourmash for proteins 
nohup python3 dnds_using_containment_index.py --ref_fasta ../../data/uniprotkb.fasta --query_fasta ../../data/500_sequences.fna --predict frame --k 7,14,21,28,35,42,49,56,63,70 --wd ../../data/ > log.txt 2>&1 &
