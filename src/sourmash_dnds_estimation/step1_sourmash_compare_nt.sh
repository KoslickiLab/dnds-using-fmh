#!/bin/bash
set -eoux pipefail

### Sketch and comapre nucleotide sequences for dN/dS estimation

#working directories for data and result output
data=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/HIT000324409/
#data=$1
sigs=${data}
wd=${data}sourmash_compare_nt/

mkdir ${wd}
#fasta files
ref=${data}ref_HIT000324409.fastn
samples=${data}queries_HIT000324409.fastn

#signature output filenames
ref_output=${sigs}ref_dna.sig 
samp_output=${sigs}queries_dna.sig.zip

#other parameters
ref_scaled=1 #scale factor for reference
samp_scaled=1 #scale factor for query

#k is k*3 for nt when comparing between protein and nt
sourmash sketch dna -p k=15,k=18,k=21,k=24,k=27,k=30,k=33,k=36,k=39,k=42,k=45,k=60,k=75,k=90,scaled=$ref_scaled $ref -o $ref_output
sourmash sketch dna -p k=15,k=18,k=21,k=24,k=27,k=30,k=33,k=36,k=39,k=42,k=45,k=60,k=75,k=90,scaled=$samp_scaled $samples --singleton -o $samp_output


for K in 15 18 21 24 27 30 33 36 39 42 45 60 75 90
do
    nohup sourmash compare $ref_output $samp_output --containment --dna --o ${wd}compare${K}.mat --csv ${wd}compare${K}.csv  --ksize $K > ${wd}compare${K}.txt 2>&1 & #assign to diff log files
    nohup sourmash compare $ref_output $samp_output --ani --dna --o ${wd}compare_ani${K}.mat --csv ${wd}compare_ani${K}.csv  --ksize $K > ${wd}compare_ani${K}.txt 2>&1 & #assign to diff log files

done


