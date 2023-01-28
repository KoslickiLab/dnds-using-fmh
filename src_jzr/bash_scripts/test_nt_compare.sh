#!/bin/bash
set -eoux pipefail

#files for a simple sequence
#wd=/data/jzr5814/data/tests/nt/simple/
#files for a single sequence
#wd=/data/jzr5814/data/tests/nt/single_seq/
#files for a 9 sequences
wd=/data/jzr5814/data/tests/nt/9_seqs/
ref=${wd}ref.fna
samples=${wd}query.fna 
ref_output=${wd}ref.sig 
samp_output=${wd}query.sig.zip
ref_scaled=1
samp_scaled=1
tbp=1

sourmash sketch dna -p k=7,k=14,k=21,k=28,k=35,k=42,k=49,k=56,k=63,k=70,scaled=$ref_scaled $ref -o $ref_output

sourmash sketch dna -p k=7,k=14,k=21,k=28,k=35,k=42,k=49,k=56,k=63,k=70,scaled=$samp_scaled $samples --singleton -o $samp_output

for K in 7 14 21 28 35 42 49 56 62 70 
do
    nohup sourmash compare $ref_output $samp_output --containment --dna --o ${wd}compare${K}.mat --ksize $K > ${wd}compare${K}.txt 2>&1 & #assign to diff log files

done


