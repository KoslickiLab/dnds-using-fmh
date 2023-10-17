#!/bin/bash
set -eoux pipefail

wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_genome_sketches'
format='cds.fna'

cd ${wd}
filelist=*.${format}

for file in $filelist; do
    #echo $file

    SUBSTRING1=$(echo $file | cut -d'.' -f 1,2)

    echo ${SUBSTRING1}

    cat ${file} | sed 's/>lcl|/>/' | sed 's/ .*//' > ${wd}/${SUBSTRING1}.cds.name_change.fna

    echo 'it worked'

    transeq -sequence ${wd}/${SUBSTRING1}.cds.name_change.fna -outseq ${wd}/${SUBSTRING1}.cds.temp.faa

    cat ${wd}/${SUBSTRING1}.cds.temp.faa | sed 's/_[0-9]*$//' > ${wd}/${SUBSTRING1}.cds.name_change.faa

done

rm -rf *temp*
