#!/bin/bash
set -eoux pipefail

### Sketch and comapre protein sequences for dN/dS estimation

#working directories for data and result output
data=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_0.1_ksizes_5_30_500_cdn_differences/
sigs=${data}
wd=${data}sourmash_compare_protein/

mkdir ${wd}

#files for a 9 sequences
ref=${data}translated_ref_seq.faa
samples=${data}translated_mutated_queries_seq.faa
#ref=${data}2499_prot_ref_seq.faa
#samples=${data}2499_prot_mutated_queries_seq.faa
#ref=${data}K03427_1452nt_ref_seq.fna
#samples=${data}K03427_1452nt_query_seqs.fna

#signature output filenames
ref_output=${sigs}ground_truth_ref_translate.sig 
samp_output=${sigs}ground_truth_mutated_queries_translate.sig.zip
#ref_output=${sigs}2499_prot_ref_seq.sig 
#samp_output=${sigs}2499_prot_mutated_queries_seq.sig.zip

#other parameters
ref_scaled=1 #scale factor for reference
samp_scaled=1 #scale factor for query

#Input of ref protein sequencess 
sourmash sketch protein -p k=5,k=6,k=7,k=8,k=9,k=10,k=15,k=20,k=25,k=30,scaled=$ref_scaled $ref -o $ref_output

#Input of mutated protein sequencess
sourmash sketch protein -p k=5,k=6,k=7,k=8,k=9,k=10,k=15,k=20,k=25,k=30,scaled=$samp_scaled $samples --singleton -o $samp_output

for K in 5 6 7 8 9 10 15 20 25 30
do

    nohup sourmash compare $ref_output $samp_output --containment --protein --o ${wd}compare${K}.mat --csv ${wd}compare${K}.csv --ksize $K > ${wd}compare${K}.txt 2>&1 & #assign to diff log files

done