#!/bin/bash

rscript=/data/jzr5814/repositories/dnds-using-fmh/src/orthologr_tool/analysis.R

wd='/data/jzr5814/orthologr_analysis/HIT000324409_pairwise'
format='fna'
mkdir ${wd}/dnds_results_via_orthologr
query=${wd}/sequences_with_no_line_breaks.fna

cd $wd/individual_fasta_files
filelist=*.${format}
for ref in $filelist; do
    echo $ref
    SUBSTRING=$(echo $ref | cut -d'.' -f 1)
    
    #nohup Rscript --vanilla $rscript ${wd}/individual_fasta_files/${ref} ${query} ${wd}/dnds_results_via_orthologr/${SUBSTRING}.csv > ${wd}/${SUBSTRING}.log 2>&1 &
done
