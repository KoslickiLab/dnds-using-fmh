#!/bin/bash
set -eoux pipefail

wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/marine_bacteria'

ref_dna=${wd}/signatures/marine_bacteria.dna.sig.gzip
query_dna=${wd}/signatures/marine_bacteria.dna.sig.gzip

ref_prot=${wd}/signatures/marine_bacteria.prot.sig.gzip
query_prot=${wd}/signatures/marine_bacteria.prot.sig.gzip

cd ${wd}
mkdir compare_protein
mkdir compare_dna

output_compare_protein=${wd}/compare_protein
output_compare_dna=${wd}/compare_dna

for ksize in 5 7 10 15 20; do
    nt_k=$(( 3*ksize ))
    
    nohup time sourmash compare $ref_dna $query_dna --containment --dna --ksize $nt_k --o ${output_compare_dna}/compare.dna.${nt_k}.mat --csv ${output_compare_dna}/compare.dna.${nt_k}.csv > ${wd}/compare.dna.${nt_k}.log 2>&1 & #assign to diff log files
    nohup time sourmash compare $ref_prot $query_prot --containment --protein --ksize $ksize --o ${output_compare_protein}/compare.prot.${ksize}.mat --csv ${output_compare_protein}/compare.prot.${ksize}.csv > ${wd}/compare.prot.${ksize}.log 2>&1 & #assign to diff log files

done
