#!/bin/bash
set -eoux pipefail

#working directories for data and result output
wd=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/HIT000324409_pairwise_take_2_with_CDS_transcripts/ #output dNdS results to this directory
nt_compare_wd=${wd}sourmash_compare_nt/ #DNA compare
prot_compare_wd=${wd}sourmash_compare_protein/ #protein/translate compare

for k in 5 6 7 8 9 10 11 12 13 14 15 20 25 30
do
declare -i nt_k=k*3

nt_containment=${nt_compare_wd}compare${nt_k}.csv
protein_containment=${prot_compare_wd}compare${k}.csv

out=dnds_constant_${k}.csv

nohup python3 dnds_CI.py --nt ${nt_containment} --protein ${protein_containment} --k ${k} --o ${out} --wd ${wd} > ${wd}dnds_CI_log.txt 2>&1 &

done

