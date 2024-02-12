#!/bin/bash
set -eoux pipefail

#Parameters
wd=/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gKaKs_paper
klist='5,7,9,11,13'
scaled=1
cds_input_list=${wd}/fasta_files.txt
out=gKaKs_paper
multiple=yes #choose yes if you have multiple fastn files

#Run FMH Omega
python3 script_fmh_omega.py --cds_input ${cds_input_list} --scaled_input ${scaled} --klist ${klist} --multiple ${multiple} --outname ${out} --working_dir ${wd}