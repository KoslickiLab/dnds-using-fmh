#!/bin/bash

wd='/data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise_take_2_with_CDS_transcripts'
format='fa'
#mkdir ${wd}/clustalw_protein_pairwise
pairwise='pairwise_cds_fna_aln'
mkdir ${wd}/${pairwise}

cd $wd/pairwise_cds_fna
filelist=*.${format}
for file1 in $filelist; do
    SUBSTRING1=$(echo $file1 | cut -d'.' -f 1)
    nohup clustalw -INFILE=${file1} -type=DNA -OUTFILE=${wd}/${pairwise}/${SUBSTRING1}.aln -OUTPUT=FASTA > ${wd}/${pairwise}/${SUBSTRING1}_aln.log 2>&1 &
done
