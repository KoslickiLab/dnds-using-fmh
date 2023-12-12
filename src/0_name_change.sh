#!/bin/bash
set -eoux pipefail

wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/nakano_sequences'
format='fna'

cd ${wd}
filelist=*.${format}
cd samples
for file in $filelist; do
    #echo $file

    SUBSTRING1=$(echo $file | cut -d'.' -f 1,2)

    echo ${SUBSTRING1}

    #cat ${file} | sed 's/>lcl|/>/' | sed 's/ .*//' > ${wd}/${SUBSTRING1}.name_change.fna

    cat ${file} | sed 's/>* /_/' | sed 's/ .*//' > ${wd}/${SUBSTRING1}.name_change.fna

    echo 'it worked'

done

rm -rf *temp*
