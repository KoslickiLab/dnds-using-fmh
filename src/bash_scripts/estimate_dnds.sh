#!/bin/bash
set -eoux pipefail

wd='/data/jzr5814/data/tests/dnds_test/'

for k in 5 10 15 20
do

nt_containment=${wd}'nt/9_seqs/compare'${k}'.csv'
protein_containment=${wd}'prot/9_seqs/compare'${k}'.csv'

out=${wd}'dNdS_'${k}.csv

nohup python3 /data/jzr5814/dnds-using-fmh/src_jzr/dnds_CI.py --nt ${nt_containment} --protein ${protein_containment} --k ${k} --o ${out} > ${wd}dnds_CI_log.txt 2>&1 &

done
