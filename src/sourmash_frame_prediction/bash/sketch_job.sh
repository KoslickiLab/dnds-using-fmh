#!/bin/bash
set -eoux pipefail

ref=/data/jzr5814/data/uniprotkb.fasta
samples=/data/jzr5814/data/frame_analysis_using_bash_script/query_frames.faa #found in /data/jzr5814/data/frame_analysis_using_bash_script
ref_scaled=100
samp_scaled=1
wd=/data/jzr5814/data/frame_analysis_using_bash_script/compare_by_multiple_jobs/test2/
ref_output=${wd}uniprotkb_test2.sig #found in /data/jzr5814/data/
samp_output=${wd}query_frames_test2.sig.zip

nohup sourmash sketch protein -p k=21,k=22,k=23,k=24,k=25,k=26,k=27,k=28,scaled=$ref_scaled $ref -o $ref_output > ${wd}sketch_ref_log.txt 2>&1 &

nohup sourmash sketch protein -p k=21,k=22,k=23,k=24,k=25,k=26,k=27,k=28,scaled=$samp_scaled $samples --singleton -o $samp_output > ${wd}sketch_protein_log.txt 2>&1 &