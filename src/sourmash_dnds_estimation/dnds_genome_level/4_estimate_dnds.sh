#!/bin/bash
set -eoux pipefail

#working directories for data and result output
wd=/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_cds_genome/ #output dNdS results to this directory

for k in 5 7 10 15 20
do
declare -i nt_k=k*3

nt_containment=${wd}/compare_dna/compare.dna.${nt_k}.csv
protein_containment=${wd}/compare_protein/compare.prot.${k}.csv

out=dnds_constant_${k}.csv

python3 dnds_CI.py --nt ${nt_containment} --protein ${protein_containment} --k ${k} --o ${out} --wd ${wd}

done

echo ${wd}dnds_constant_*

awk 'NR == 1 || FNR > 1' ${wd}dnds_constant_5.csv ${wd}dnds_constant_7.csv ${wd}dnds_constant_10.csv ${wd}dnds_constant_15.csv ${wd}dnds_constant_20.csv > ${wd}all_dnds_constant.csv
