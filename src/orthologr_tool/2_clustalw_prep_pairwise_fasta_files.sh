#!/bin/bash

wd='/data/jzr5814/orthologr_analysis/HIT000324409_pairwise'
format='cds'
mkdir ${wd}/prep_pairwise_${format}
mkdir ${wd}/prep_pairwise_faa_from_cds


cd $wd/individual_${format}_files
#filelist=*.${format}
filelist=*.fa
for file1 in $filelist; do
    for file2 in $filelist; do
        fasta_file1=${file1}
        fasta_file2=${file2}

        SUBSTRING1=$(echo $file1 | cut -d'.' -f 1)
        SUBSTRING2=$(echo $file2 | cut -d'.' -f 1)

        if [ ${fasta_file1} != ${fasta_file2} ]; then
            cat ${wd}/individual_${format}_files/${fasta_file1} ${wd}/individual_${format}_files/${fasta_file2} > ${wd}/prep_pairwise_${format}/${SUBSTRING1}_${SUBSTRING2}.${format} #save in a different directory
            transeq -sformat pearson ${wd}/prep_pairwise_${format}/${SUBSTRING1}_${SUBSTRING2}.${format} ${wd}/prep_pairwise_faa_from_cds/${SUBSTRING1}_${SUBSTRING2}.faa 
        fi
    done
done

