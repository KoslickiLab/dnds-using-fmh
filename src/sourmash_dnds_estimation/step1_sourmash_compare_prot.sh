#!/bin/bash
set -eoux pipefail

### Sketch and comapre protein sequences for dN/dS estimation

#working directories for data and result output
data=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/
sigs=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/ #sketch sig output directory
wd=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/sourmash_compare_protein/

mkdir ${wd}

#files for a 9 sequences
#ref=${data}ref.fna
#samples=${data}query.fna
#ref=${data}ground_truth_ref_used_for_dNdS.fna
#samples=${data}ground_truth_mutated_queries_used_for_dNdS.fna
ref=${data}10000_prot_ref_seq.faa
samples=${data}10000_prot_mutated_queries_seq.faa

#signature output filenames
#ref_output=${sigs}ref.sig 
#samp_output=${sigs}query.sig.zip
#ref_output=${sigs}ground_truth_ref_translate.sig 
#samp_output=${sigs}ground_truth_mutated_queries_translate.sig.zip
ref_output=${sigs}10000_prot_ref_seq.sig 
samp_output=${sigs}10000_prot_mutated_queries_seq.sig.zip

#other parameters
ref_scaled=1 #scale factor for reference
samp_scaled=1 #scale factor for query

#Input of ref protein sequencess 
sourmash sketch protein -p k=7,scaled=$ref_scaled $ref -o $ref_output

#Input of mutated protein sequencess
sourmash sketch protein -p k=7,scaled=$samp_scaled $samples --singleton -o $samp_output

#for K in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 20
for K in 7
do

    nohup sourmash compare $ref_output $samp_output --containment --protein --o ${wd}compare${K}.mat --csv ${wd}compare${K}.csv --ksize $K > ${wd}compare${K}.txt 2>&1 & #assign to diff log files

done