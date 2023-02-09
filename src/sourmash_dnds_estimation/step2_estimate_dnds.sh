#!/bin/bash
set -eoux pipefail

#working directories for data and result output
nt_compare_wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/sourmash_compare/sourmash_compare_nt/' #DNA compare
prot_compare_wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/sourmash_compare/sourmash_compare_protein/' #protein/translate compare
wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/DNDS/' #output dNdS results to this directory

for k in 5 10 15 20
do

nt_containment=${nt_compare_wd}'compare'${k}'.csv'
protein_containment=${prot_compare_wd}'compare'${k}'.csv'

out=${wd}'dNdS_'${k}.csv

nohup python3 dnds_CI.py --nt ${nt_containment} --protein ${protein_containment} --k ${k} --o ${out} > ${wd}dnds_CI_log.txt 2>&1 &

done
