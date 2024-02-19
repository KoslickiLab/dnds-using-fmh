#!/bin/bash
set -eoux pipefail

#Parameters
wd=/data/jzr5814/sourmash_dnds_estimation/tests/test/branchwater_test/positive_subsets_0.01
klist='5,7'
scaled=1
cds_input_list=${wd}/fasta_files.txt
m="branchwater" #branchwater mode, single fasta file use sourmash sketch, multiple fasta files to use sourmash sketch
outname=gKaKs_paper

#Run FMH Omega
python3 script_fmh_omega.py --fasta_input_list ${cds_input_list} --scaled_input ${scaled} --klist ${klist} --mode ${m} --outname ${out} --working_dir ${wd}
