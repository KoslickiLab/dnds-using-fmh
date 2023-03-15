#!/bin/bash
set -eoux pipefail

### Sketch and comapre nucleotide sequences for dN/dS estimation

#working directories for data and result output
data=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/
sigs=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/
wd=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/sourmash_compare_nt/

mkdir ${wd}
#files for a 9 sequences
#ref=${data}ref.fna
#samples=${data}query.fna
#ref=${data}ground_truth_ref_used_for_dNdS.fna
#samples=${data}ground_truth_mutated_queries_used_for_dNdS.fna
ref=${data}10000_nt_ref_seq.fna
samples=${data}10000_nt_mutated_queries_seq.fna

#signature output filenames
#ref_output=${sigs}ref.sig 
#samp_output=${sigs}query.sig.zip
#ref_output=${sigs}ground_truth_ref_dna.sig 
#samp_output=${sigs}ground_truth_mutated_queries_dna.sig.zip
ref_output=${sigs}10000_nt_ref_seq.sig 
samp_output=${sigs}10000_nt_ref_seq.sig.zip

#other parameters
ref_scaled=1 #scale factor for reference
samp_scaled=1 #scale factor for query

#k is k*3 for nt when comparing between protein and nt
sourmash sketch dna -p k=21,scaled=$ref_scaled $ref -o $ref_output
#sourmash sketch dna -p k=2,k=3,k=4,k=5,k=6,k=7,k=8,k=9,k=10,k=11,k=12,k=13,k=14,k=15,k=20,scaled=$ref_scaled $ref -o $ref_output

sourmash sketch dna -p k=21,scaled=$samp_scaled $samples --singleton -o $samp_output
#sourmash sketch dna -p k=2,k=3,k=4,k=5,k=6,k=7,k=8,k=9,k=10,k=11,k=12,k=13,k=14,k=15,k=20,scaled=$samp_scaled $samples --singleton -o $samp_output


#for K in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 20 
for K in 21
do
    nohup sourmash compare $ref_output $samp_output --containment --dna --o ${wd}compare${K}.mat --csv ${wd}compare${K}.csv  --ksize $K > ${wd}compare${K}.txt 2>&1 & #assign to diff log files

done


