#!/bin/bash

wd='/data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise_take_2_with_CDS_transcripts'
format='fa'
#mkdir ${wd}/clustalw_protein_pairwise
pairwise='pairwise_cds_fna'
mkdir ${wd}/${pairwise}

cd $wd/individual_cds_fna
filelist=*.${format}
for file1 in $filelist; do
    for file2 in $filelist; do
        SUBSTRING1=$(echo $file1 | cut -d'.' -f 1)
        SUBSTRING2=$(echo $file2 | cut -d'.' -f 1)
        if [${SUBSTRING1} != ${SUBSTRING2}]; then
            cat ${wd}/individual_cds_fna/$file1 ${wd}/individual_cds_fna/${file2} > ${wd}/${pairwise}/${SUBSTRING1}_${SUBSTRING2}.fa #save in a different directory
        fi
    done
done

#rm ${wd}/individual_fasta_files/*
#rmdir ${wd}/individual_fasta_files
