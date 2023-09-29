#!/bin/bash

python kaks_calc_prep.py
python clustalw_job_prep.py

wd='/data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise'

mkdir ${wd}/clustalw_protein_pairwise

cd $wd/individual_fasta_files
filelist=*.faa
for file1 in $filelist; do
    for file2 in $filelist; do
        cat ${wd}/individual_fasta_files/$file1 ${wd}/individual_fasta_files/${file2} > ${wd}/clustalw_protein_pairwise/${file1%..faa}_${file2} #save in a different directory
    done
done

rm ${wd}/individual_fasta_files/*
rmdir ${wd}/individual_fasta_files/