#!/bin/bash
set -eoux pipefail

#files for a simple sequence
#wd=/data/jzr5814/data/tests/nt/simple/
#files for a single sequence
wd=/data/jzr5814/data/tests/nt/single_seq/
#files for a 9 sequences
#wd=/data/jzr5814/data/tests/nt/9_seqs/
ref=${wd}ref.fna
samples=${wd}query.fna 
ref_output=${wd}ref.sig 
samp_output=${wd}query.sig.zip
ref_scaled=1
samp_scaled=1
K=7

sourmash sketch dna -p k=5,k=10,k=15,k=20,scaled=$ref_scaled $ref -o $ref_output

sourmash sketch dna -p k=5,k=10,k=15,k=20,scaled=$samp_scaled $samples --singleton -o $samp_output

#nohup sourmash compare $ref_output $samp_output --containment --dna --o ${wd}compare${K}.mat --csv ${wd}compare${K}.csv  --ksize $K > ${wd}compare${K}.txt 2>&1 &

#for K in 7 14 21 28 35 42 49 56 62 70
for K in 5 10 15 20 
do
    nohup sourmash compare $ref_output $samp_output --containment --dna --o ${wd}compare${K}.mat --csv ${wd}compare${K}.csv  --ksize $K > ${wd}compare${K}.txt 2>&1 & #assign to diff log files

done


