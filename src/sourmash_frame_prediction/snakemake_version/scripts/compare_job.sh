#!/bin/bash
set -eoux pipefail

ref_output=$1
samp_output=$2
ksize=$3
wd=$4

echo $ksize
echo $wd
#for K in 7 14 21 28 35 42 49 56 63 70
for K in $ksize
do

    echo $K
    #nohup sourmash compare $ref_output $samp_output --containment --protein --o ${wd}compare${K}.mat --csv ${wd}compare${K}.csv --ksize $K > ${wd}compare${K}.txt 2>&1 & #assign to diff log files

done