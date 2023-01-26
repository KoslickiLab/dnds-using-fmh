#!/bin/bash
set -eoux pipefail

#Test run on small dataset
#nohup python3 dnds_using_containment_index.py --ref_fasta ../../data/tests/50_sequences.fna --query_fasta ../../data/tests/5_nt_sequences.fna --predict frame --k 7,14,21,28 --wd ../../data/tests/ > log.txt 2>&1 &

#Running sourmash for dna
nohup python3 dnds_using_containment_index.py --ref_fasta /data/jzr5814/data/tests/nt/simple/ref.fna --query_fasta /data/jzr5814/data/tests/nt/simple/query.fna --k 3,4,5,6 --wd /data/jzr5814/data/tests/nt/simple/ --moltype dna --translate no > nt_log.txt 2>&1 &

#Running sourmash for protein
nohup python3 dnds_using_containment_index.py --ref_fasta /data/jzr5814/data/tests/prot/simple/ref.fna --query_fasta /data/jzr5814/data/tests/prot/simple/query.fna --predict frame --k 3,4,5,6 --wd /data/jzr5814/data/tests/prot/simple/ --moltype protein --translate yes > prot_log.txt 2>&1 &
