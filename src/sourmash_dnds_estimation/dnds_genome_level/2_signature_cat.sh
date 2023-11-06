#!/bin/bash
time -v
set -eoux pipefail

start=`date +%s`

wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_nt_genome_redo'

cd ${wd}
mkdir signature_cat

nohup sourmash signature cat ${wd}/signatures/*.dna.sig.gzip -o ${wd}/signature_cat/dna.cat.sig.gzip > ${wd}/cat_dna.log 2>&1 &
nohup sourmash signature cat ${wd}/signatures/*.prot.sig.gzip -o ${wd}/signature_cat/prot.cat.sig.gzip > ${wd}/cat_protein.log 2>&1 &

