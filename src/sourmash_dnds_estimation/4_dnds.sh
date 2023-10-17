#!/bin/bash
set -eoux pipefail

#working directories for data and result output
wd=/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_genome_sketches/test_small/ #output dNdS results to this directory

for k in 5 7 10 15 20
do
declare -i nt_k=k*3

nt_containment=${wd}compare.dna.${nt_k}.csv
protein_containment=${wd}compare.protein.${k}.csv

out=dnds_constant_${k}.csv

nohup python3 dnds_CI.py --nt ${nt_containment} --protein ${protein_containment} --k ${k} --o ${out} --wd ${wd} > ${wd}dnds_CI_log.txt 2>&1 &

done

echo ${wd}dnds_constant_*

awk 'NR == 1 || FNR > 1' ${wd}dnds_constant_* > ${wd}dnds_constant_all.csv
