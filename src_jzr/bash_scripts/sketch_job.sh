#!/bin/bash
set -eoux pipefail

ref=../../data/uniprotkb.fasta
samples=../../data/query_frames.faa #found in /data/jzr5814/data/frame_analysis_using_bash_script
ref_scaled=100
samp_scaled=1
wd=/data/jzr5814/data/frame_analysis_using_bash_script/compare_by_multiple_jobs/test2/
ref_output=${wd}uniprotkb_test2.sig #found in /data/jzr5814/data/
samp_output=${wd}query_frames_test2.sig.zip

nohup sourmash sketch protein -p k=7,k=14,k=15,k=16,k=17,k=18,k=19,k=20,k=21,scaled=$ref_scaled $ref -o $ref_output > ${wd}sketch_ref_log.txt 2>&1 &

nohup sourmash sketch protein -p k=7,k=14,k=15,k=16,k=17,k=18,k=19,k=20,k=21,scaled=$samp_scaled $samples --singleton -o $samp_output > ${wd}sketch_protein_log.txt 2>&1 &
