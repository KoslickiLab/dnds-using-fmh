#!/bin/bash
set -eoux pipefail

#s=1
#for ksize in 7 14 21 28 35 42 49 56 63 70
#do
#echo index protein $ksize
#sourmash index /data/jzr5814/data/frame_analysis_using_bash_script/querydb_${ksize} /data/jzr5814/data/frame_analysis_using_bash_script/query_frames_sigs --protein --k ${ksize} --scaled ${s}
#done

wd=/data/jzr5814/data/frame_analysis_using_bash_script/search/test2/
ref_output=/data/jzr5814/data/uniprotkb.sig
t=0.0 # minimum threshold for matching
n=0 #minimum number of results reported

for ksize in 7 14 21 28 35 #42 49 56 63 70
do
echo search protein $ksize
db=/data/jzr5814/data/frame_analysis_using_bash_script/querydb_${ksize}.sbt.zip
nohup sourmash search $ref_output $db --protein --k ${ksize} --containment -t ${t} -n ${n} -o ${wd}search_${ksize}.csv > ${wd}search_log_${ksize}.txt 2>&1 & #assign to diff log files
done
