#!/bin/bash
set -eoux pipefail

### Sketch and comapre nucleotide sequences for dN/dS estimation

#working directories for data and result output
data=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_0.1_ksizes_5_10/
sigs=${data}
wd=${data}sourmash_compare_nt/

mkdir ${wd}
#fasta files
ref=${data}ground_truth_ref_used_for_dNdS.fna
samples=${data}ground_truth_mutated_queries_used_for_dNdS.fna
#ref=${data}10000_nt_ref_seq.fna
#samples=${data}10000_nt_mutated_queries_seq.fna

#signature output filenames
ref_output=${sigs}ground_truth_ref_dna.sig 
samp_output=${sigs}ground_truth_mutated_queries_dna.sig.zip
#ref_output=${sigs}10000_nt_ref_seq.sig 
#samp_output=${sigs}10000_nt_ref_seq.sig.zip

#other parameters
ref_scaled=1 #scale factor for reference
samp_scaled=1 #scale factor for query

#k is k*3 for nt when comparing between protein and nt
#sourmash sketch dna -p k=21,scaled=$ref_scaled $ref -o $ref_output
sourmash sketch dna -p k=30,k=33,k=36,k=39,k=42,k=45,scaled=$ref_scaled $ref -o $ref_output

#sourmash sketch dna -p k=21,scaled=$samp_scaled $samples --singleton -o $samp_output
sourmash sketch dna -p k=30,k=33,k=36,k=39,k=42,k=45,scaled=$samp_scaled $samples --singleton -o $samp_output


for K in 30 33 36 39 42 45
do
    nohup sourmash compare $ref_output $samp_output --containment --dna --o ${wd}compare${K}.mat --csv ${wd}compare${K}.csv  --ksize $K > ${wd}compare${K}.txt 2>&1 & #assign to diff log files

done


