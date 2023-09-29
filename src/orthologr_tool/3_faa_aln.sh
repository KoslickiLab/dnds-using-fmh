#!/bin/bash

#set wroking directory to find files
wd='/data/jzr5814/orthologr_analysis/HIT000324409_pairwise'

#create new folder to store alignment files
mkdir ${wd}/pairwise_faa_aln_new

#go into faa files prepared for aln
cd $wd/prep_pairwise_faa_from_cds

#get list of files in current directory (what is the current directory?)
filelist=*.faa

#for loop through these files
for input_file in $filelist; do

    SUBSTRING=$(echo $input_file | cut -d'.' -f 1)
    output_file=${wd}/pairwise_faa_aln_new/${SUBSTRING}.aln
    log_file=${wd}/pairwise_faa_aln_new/${SUBSTRING}.log
    
    #align using clustalw
    nohup clustalw -INFILE=$input_file -type=PROTEIN -OUTFILE=$output_file -OUTPUT=FASTA > $log_file 2>&1 &
    done