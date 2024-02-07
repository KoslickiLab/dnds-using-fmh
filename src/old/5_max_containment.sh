#!/bin/bash
set -eoux pipefail

#working directories for data and result output
wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/nakano_sequences/'

#get the max containment
for k in 5 7 10 15 20
do
declare -i nt_k=k*3

dnds_constant=${wd}dnds_constant_${k}.csv

out=dnds_max_nt_containments_${k}_with_selection_pressure.csv
python3 max_nt_containment.py --dnds_constant ${dnds_constant} --k ${k} --o ${out} --wd ${wd}

out=dnds_max_protein_containments_${k}_with_selection_pressure.csv
python3 max_protein_containment.py --dnds_constant ${dnds_constant} --k ${k} --o ${out} --wd ${wd}

done

#concatenate all max containments of kmers into one file
awk 'NR == 1 || FNR > 1' ${wd}dnds_max_nt_containments_5_with_selection_pressure.csv ${wd}dnds_max_nt_containments_7_with_selection_pressure.csv ${wd}dnds_max_nt_containments_10_with_selection_pressure.csv ${wd}dnds_max_nt_containments_15_with_selection_pressure.csv ${wd}dnds_max_nt_containments_20_with_selection_pressure.csv > ${wd}all_dnds_max_nt_containments_with_selection_pressure.csv
awk 'NR == 1 || FNR > 1' ${wd}dnds_max_protein_containments_5_with_selection_pressure.csv ${wd}dnds_max_protein_containments_7_with_selection_pressure.csv ${wd}dnds_max_protein_containments_10_with_selection_pressure.csv ${wd}dnds_max_protein_containments_15_with_selection_pressure.csv ${wd}dnds_max_protein_containments_20_with_selection_pressure.csv > ${wd}all_dnds_max_protein_containments_with_selection_pressure.csv

#fix column names for figure generation
max_dnds_constant=${wd}all_dnds_max_nt_containments_with_selection_pressure.csv
out=all_dnds_max_nt_containments_with_selection_pressure_fixed_columns_names.csv
python3 fix_column_names_for_merging.py --max_dnds_constant ${max_dnds_constant} --o ${out} --wd ${wd}

max_dnds_constant=${wd}all_dnds_max_protein_containments_with_selection_pressure.csv
out=all_dnds_max_protein_containments_with_selection_pressure_fixed_columns_names.csv
python3 fix_column_names_for_merging.py --max_dnds_constant ${max_dnds_constant} --o ${out} --wd ${wd}
