#!/bin/bash
set -eoux pipefail

#Parameters
wd=/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gtdb_protein_rep_pairwise_archaea
k='5'
scaled=1
translate_cds='no' #indicate with yes or no
m="bwpair" #branchwater mode, single fasta file use sourmash sketch, multiple fasta files to use sourmash sketch

#Note that if using branchwater mode, makesure working directory is at the beginning of filename. THIS IS NOT REQUIERED FOR OTHER MODES
cds_input_list=${wd}/dataset.csv
out=demo_branchwater

#Run FMH Omega
python3 script_fmh_omega.py --fasta_input_list ${cds_input_list} --scaled_input ${scaled} --ksize ${k} --mode ${m} --translate ${translate_cds} --outname ${out} --directory ${wd} 
