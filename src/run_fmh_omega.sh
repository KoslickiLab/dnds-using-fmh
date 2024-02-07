#!/bin/bash
set -eoux pipefail

#Parameters
wd=/data/jzr5814/sourmash_dnds_estimation/tests/test/testing_new_code_on_create_sequence_using_NG_assumption_0.01/positive
klist='5,7'
scaled=1
cds_input_file=${wd}/positive_selection_queries_10002_0.01_small.fna
out=testing_new_code
singles=yes

#Run FMH Omega
python3 fmh_omega_from_cds.py --cds_input ${cds_input_file} --scaled_input ${scaled} --klist ${klist} --singleton ${singles} --outname ${out} --working_dir ${wd}