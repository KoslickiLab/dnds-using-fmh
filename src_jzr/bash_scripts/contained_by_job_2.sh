#!/bin/bash
set -eoux pipefail

for ksize in 7 14 21 28 35 42 49 56 63 70 
do
nohup python3 /data/jzr5814/dnds-using-fmh/src_jzr/sourmash_api_frame_cont_2.py --ksize ${ksize} > /data/jzr5814/data/frame_analysis_using_bash_script/contained_by_multiple_jobs/contained_by_${ksize}_log.txt 2>&1 &
done
