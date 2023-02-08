#!/bin/bash
set -eoux pipefail

#results from this are found in /data/jzr5814/data/frame_analysis_using_bash_script
#Decided to run abashscript because the python program has issues

ref=../../data/uniprotkb.fasta
samples=../../data/query_frames.faa #found in /data/jzr5814/data/frame_analysis_using_bash_script
ref_output=/data/jzr5814/data/uniprotkb.sig #found in /data/jzr5814/data/
samp_output=/data/jzr5814/data/frame_analysis_using_bash_script/query_frames.sig.zip
ref_scaled=100
samp_scaled=1
tbp=0 #threshold of 0 should report all matched and unmatched hashes
wd=/data/jzr5814/data/frame_analysis_using_bash_script/

#sourmash sketch protein -p k=7,k=14,k=21,k=28,k=35,k=42,k=49,k=56,k=63,k=70,scaled=$ref_scaled $ref -o $ref_output

#sourmash sketch protein -p k=7,k=14,k=21,k=28,k=35,k=42,k=49,k=56,k=63,k=70,scaled=$samp_scaled $samples --singleton -o $samp_output

for K in 7 14 21 28 35 42 49 56 63 70
do
	echo $K
	nohup sourmash prefetch $ref_output $samp_output --protein --o ${wd}prefetch_res_${K}.csv --threshold-bp $tbp --ksize $K > ${wd}test_prefetch${K}.txt 2>&1 & #assign to diff log files
done

python3 ../frame_analysis.py --wd ${wd}

