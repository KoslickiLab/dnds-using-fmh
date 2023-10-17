#!/bin/bash
set -eoux pipefail

### Sketch and comapre protein sequences for dN/dS estimation

#working directories for data and result output
data=/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_with_singleton/
sigs=${data}
wd=${data}sourmash_compare_protein/

mkdir ${wd}

#files for a 9 sequences
ref=${data}ref_ecoli_10strains_name_change.faa
samples=${data}ecoli_10strains_name_change.faa

#signature output filenames
ref_output=${sigs}ref_translate.sig 
samp_output=${sigs}queries_translate.sig.zip

#other parameters
ref_scaled=1 #scale factor for reference
samp_scaled=1 #scale factor for query

#Input of ref protein sequencess 
sourmash sketch protein -p k=5,k=6,k=7,k=8,k=9,k=10,k=11,k=12,k=13,k=14,k=15,k=20,k=25,k=30,scaled=$ref_scaled $ref --singleton -o $ref_output
#Input of mutated protein sequencess
sourmash sketch protein -p k=5,k=6,k=7,k=8,k=9,k=10,k=11,k=12,k=13,k=14,k=15,k=20,k=25,k=30,scaled=$samp_scaled $samples --singleton -o $samp_output

for K in 5 6 7 8 9 10 11 12 13 14 15 20 25 30
do

    nohup sourmash compare $ref_output $samp_output --containment --protein --o ${wd}compare${K}.mat --csv ${wd}compare${K}.csv --ksize $K > ${wd}compare${K}.txt 2>&1 & #assign to diff log files
    nohup sourmash compare $ref_output $samp_output --ani --protein --o ${wd}compare_ani${K}.mat --csv ${wd}compare_ani${K}.csv  --ksize $K > ${wd}compare_ani${K}.txt 2>&1 & #assign to diff log files

done