#!/bin/bash
set -eoux pipefail

#files for a simple sequence
#wd=/data/jzr5814/data/tests/prot/simple/
#files for a single sequence
#wd=/data/jzr5814/data/tests/prot/single_seq/
#files for a 9 sequences
wd=/data/jzr5814/data/tests/prot/9_seqs/
ref=${wd}ref.fna
samples=${wd}query.fna 
ref_output=${wd}ref.sig 
samp_output=${wd}query.sig.zip
ref_scaled=1
samp_scaled=1
tbp=1

sourmash sketch protein -p k=7,scaled=$ref_scaled $ref -o $ref_output

sourmash sketch protein -p k=7,scaled=$samp_scaled $samples --singleton -o $samp_output

K=7
nohup sourmash compare $ref_output $samp_output --containment --protein --o ${wd}compare${K}.mat --ksize $K > ${wd}compare${K}.txt 2>&1 & #assign to diff log files

