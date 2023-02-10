#!/bin/bash
set -eoux pipefail

#files for a simple sequence
wd=/data/jzr5814/data/frame_analysis_using_bash_script/compare_by_multiple_jobs/
ref_output=/data/jzr5814/data/uniprotkb.sig
samp_output=/data/jzr5814/data/frame_analysis_using_bash_script/query_frames.sig.zip
ref_scaled=100
samp_scaled=1

for K in 7 14 21 28 35 42 49 56 63 70 
do

    nohup sourmash compare $ref_output $samp_output --containment --protein --o ${wd}compare${K}.mat --csv ${wd}compare${K}.csv --ksize $K > ${wd}compare${K}.txt 2>&1 & #assign to diff log files

done