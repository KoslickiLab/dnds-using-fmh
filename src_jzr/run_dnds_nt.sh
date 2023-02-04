#!/bin/bash
set -eoux pipefail

ref=ref.fna
samples=query.fna
ref_scaled=100
samp_scaled=1
tbp=1
wd_dna=/data/jzr5814/data/tests/nt/simple/
wd_prot=/data/jzr5814/data/tests/prot/simple/

#Test run on small dataset
#nohup python3 dnds_using_containment_index.py --ref_fasta ../../data/tests/50_sequences.fna --query_fasta ../../data/tests/5_nt_sequences.fna --predict frame --k 7,14,21,28 --wd ../../data/tests/ > log.txt 2>&1 &

#Running sourmash for dna
nohup python3 dnds_using_containment_index.py --ref_fasta ${wd_dna}ref.fna --query_fasta ${wd_dna}query.fna --k 3,4,5,6 --wd ${wd_dna} --analyze dna > nt_log.txt 2>&1 &

#Running sourmash for protein
nohup python3 dnds_using_containment_index.py --ref_fasta ${wd_prot}ref.fna --query_fasta ${wd_prot}query.fna --predict frame --k 3,4,5,6 --wd ${wd_prot} --analyze protein > prot_log.txt 2>&1 &

#Running sourmash for dna
#nohup python3 dnds_using_containment_index.py --ref_fasta /data/jzr5814/data/tests/nt/single/ref.fna --query_fasta /data/jzr5814/data/tests/nt/single/query.fna --k 3,4,5,6 --wd /data/jzr5814/data/tests/nt/single/ --moltype dna --translate no > nt_single_log.txt 2>&1 &

#Running sourmash for protein
#nohup python3 dnds_using_containment_index.py --ref_fasta /data/jzr5814/data/tests/prot/single/ref.fna --query_fasta /data/jzr5814/data/tests/prot/single/query.fna --predict frame --k 3,4,5,6 --wd /data/jzr5814/data/tests/prot/single/ --moltype protein --translate yes > prot_single_log.txt 2>&1 &

