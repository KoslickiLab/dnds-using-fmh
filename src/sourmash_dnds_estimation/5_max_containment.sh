#!/bin/bash
set -eoux pipefail

#working directories for data and result output
wd=/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/HIT000324409_pairwise_take_2_with_CDS_transcripts_change_names/ #output dNdS results to this directory

for k in 5 7 10 15 20
do
declare -i nt_k=k*3

dnds_constant=${wd}dnds_constant_${k}.csv

out=dnds_max_containments_${k}_with_selection_pressure.csv

python3 maxcontainment.py --dnds_constant ${dnds_constant} --k ${k} --o ${out} --wd ${wd}

done

awk 'NR == 1 || FNR > 1' ${wd}dnds_max_containments_5_with_selection_pressure.csv ${wd}dnds_max_containments_7_with_selection_pressure.csv ${wd}dnds_max_containments_10_with_selection_pressure.csv ${wd}dnds_max_containments_15_with_selection_pressure.csv ${wd}dnds_max_containments_20_with_selection_pressure.csv > ${wd}all_dnds_max_containments_with_selection_pressure.csv
