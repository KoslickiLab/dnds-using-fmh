#!/bin/bash

rscript=/data/jzr5814/repositories/dnds-using-fmh/src/orthologr_tool/produce_cdn_aln.R

wd='/data/jzr5814/orthologr_analysis/HIT000324409_pairwise'

mkdir ${wd}/pairwise_codon_aln_new

#query=${wd}/sequences_with_no_line_breaks.fna

cd $wd/pairwise_faa_aln_new
filelist=*.aln

for faa_aln in $filelist; do
    SUBSTRING=$(echo $faa_aln | cut -d'.' -f 1)
    #nohup Rscript --vanilla $rscript ${faa_aln} ${wd}/prep_pairwise_fna/${SUBSTRING}.fna ${wd}/pairwise_codon_aln/${SUBSTRING}.aln > ${wd}/pairwise_codon_aln/${SUBSTRING}.log 2>&1 &
    nohup Rscript --vanilla $rscript ${faa_aln} ${wd}/prep_pairwise_cds/${SUBSTRING}.cds ${wd}/pairwise_codon_aln_new/${SUBSTRING}.aln > ${wd}/pairwise_codon_aln_new/${SUBSTRING}.log 2>&1 &
done
