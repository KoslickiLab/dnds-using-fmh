#!/bin/bash
set -eoux pipefail

#Parameters
wd=/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gKaKs_paper
klist='5,7'
scaled=1
cds_input_list=${wd}/fasta_files.txt
out=gKaKs_paper
singles=no

#Run FMH Omega
python3 fmh_omega_from_cds.py --cds_input ${cds_input_list} --scaled_input ${scaled} --klist ${klist} --singleton ${singles} --outname ${out} --working_dir ${wd}