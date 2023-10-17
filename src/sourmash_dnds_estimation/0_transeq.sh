#!/bin/bash
set -eoux pipefail

wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_genome_sketches/test_small'
format='cds.fna'

cd ${wd}
filelist=*.${format}

for file in $filelist; do
    #echo $file
    SUBSTRING1=$(echo $file | cut -d'.' -f 1,2)

    cat ${file} | sed 's/>lcl|/>/' | sed 's/ .*//' > ${wd}/${SUBSTRING1}.name_change.cds.fna

    transeq ${wd}/${SUBSTRING1}.name_change.cds.fna ${wd}/${SUBSTRING1}.temp.cds.faa

    cat ${wd}/${SUBSTRING1}.temp.cds.faa | sed 's/_[0-9]*$//' > ${wd}/${SUBSTRING1}.name_change.cds.faa

done

rm -rf *temp*
