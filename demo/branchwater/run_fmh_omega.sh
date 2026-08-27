#!/bin/bash
set -eoux pipefail

#Parameters
wd=/data/jzr5814/repositories/dnds-using-fmh/src/demo/branchwater
klist='5,7'
scaled=1

#Note that if using branchwater mode, makesure working directory is at the beginning of filename. THIS IS NOT REQUIERED FOR OTHER MODES
cds_input_list=${wd}/datasets.csv

m="branchwater" #branchwater mode, single fasta file use sourmash sketch, multiple fasta files to use sourmash sketch
out=demo_branchwater

#Run FMH Omega
python3 script_fmh_omega.py --fasta_input_list ${cds_input_list} --scaled_input ${scaled} --klist ${klist} --mode ${m} --outname ${out} --directory ${wd} 
