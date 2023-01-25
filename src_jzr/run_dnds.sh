#!/bin/bash
set -eoux pipefail

nohup python3 dnds_using_containment_index.py --ref_fasta ../../data/tests/50_sequences.fna --query_fasta ../../data/tests/5_nt_sequences.fna --predict frames -k 7,14,21,28 --wd ../../data/test2/ > log.txt 2>&1 &
