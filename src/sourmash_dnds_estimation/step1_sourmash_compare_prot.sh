#!/bin/bash
set -eoux pipefail

### Sketch and comapre protein sequences for dN/dS estimation

#working directories for data and result output
data=/data/jzr5814/sourmash_dnds_estimation/tests/data/ #.faa and .fna files are found
sigs=/data/jzr5814/sourmash_dnds_estimation/tests/data/signatures/prot/ #sketch sig output directory
wd=/data/jzr5814/sourmash_dnds_estimation/tests/results/sourmash_compare/sourmash_compare_protein/ #output for sourmash comapre directory

#files for a 9 sequences
ref=${data}ref.fna
samples=${data}query.fna

#signature output filenames
ref_output=${sigs}ref.sig 
samp_output=${sigs}query.sig.zip

#other parameters
ref_scaled=1 #scale factor for reference
samp_scaled=1 #scale factor for query

#Input are protein sequencess
#sourmash sketch protein -p k=7,k=14,k=21,k=28,k=35,k=42,k=49,k=56,k=63,k=70,scaled=$ref_scaled $ref -o $ref_output
#sourmash sketch protein -p k=7,k=14,k=21,k=28,k=35,k=42,k=49,k=56,k=63,k=70,scaled=$samp_scaled $samples --singleton -o $samp_output

#Input are DNA sequences
sourmash sketch translate -p k=5,k=10,k=15,k=20,scaled=$ref_scaled $ref -o $ref_output
sourmash sketch translate -p k=5,k=10,k=15,k=20,scaled=$samp_scaled $samples --singleton -o $samp_output

for K in 5 10 15 20 
do

    nohup sourmash compare $ref_output $samp_output --containment --protein --o ${wd}compare${K}.mat --csv ${wd}compare${K}.csv --ksize $K > ${wd}compare${K}.txt 2>&1 & #assign to diff log files

done