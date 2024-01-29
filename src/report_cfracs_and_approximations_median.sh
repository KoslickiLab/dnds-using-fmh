#!/bin/bash
set -eoux pipefail

#working directories for data and result output
wd='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/fmh_dnds_sketch_translate/negative_selection_k5-15/'

for k in 5 7 9 11 13 15 21 27 33 39 45; 
do
declare -i nt_k=k*3

#fixing compare results to obtain median cfracs
nt_containment=${wd}/compare_dna/compare.dna.${k}.csv
protein_containment=${wd}/compare_protein/compare.prot.${k}.csv

python3 report_Cfracs_by_ksize.py --nt ${nt_containment} --protein ${protein_containment} --k ${k} --wd ${wd}

done

# for both nt and protein containments, make sure to run: awk 'NR == 1 || FNR > 1' nt_containment* > nt_all_cfrac.csv
awk 'NR == 1 || FNR > 1' ${wd}/ANI_approximations* > ${wd}/ANI_medians.csv
awk 'NR == 1 || FNR > 1' ${wd}/AAI_approximations* > ${wd}/AAI_medians.csv
klist='5,7,9,11,13,15,21,27,33,39,45'
python3 report_all_median_Cfracs.py --nt ${wd}/ANI_medians.csv --protein ${wd}/AAI_medians.csv --k ${klist} --wd ${wd} > log
