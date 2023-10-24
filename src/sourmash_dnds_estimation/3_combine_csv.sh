#!/bin/bash
set -eoux pipefail

wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_genome_sketches'

#delete empty csv files
find ${wd}/compare_dna/. -type f -empty -print -delete
find ${wd}/compare_protein/. -type f -empty -print -delete

for k in 5 7 10 15 20; do
    echo ${k}
    k_nt=$(( 3*k ))
    echo ${k_nt}
    #combine proteina nd dna csv files
    awk 'NR == 1 || FNR > 1' ${wd}/compare_dna/compare.dna.${k_nt}*csv | cut -d ',' -f 1,3,5 | sed 's/filename/B/' | sed 's/query_filename/A/' | awk -F ',' 'BEGIN{FS=OFS='${k}'}{print value OFS $0}' | sed "s/^${k}/${k},/" | sed "1 s/${k}/ksize/" | sed "s/similarity/containment/" | sed "s/.dna.sig.gzip//" | sed "s/.cds.name_change.fna//" > ${wd}/compare.dna.${k_nt}.csv
    awk 'NR == 1 || FNR > 1' ${wd}/compare_protein/compare.prot.${k}*csv | cut -d ',' -f 1,3,5 | sed 's/filename/B/' | sed 's/query_filename/A/' | awk -F ',' 'BEGIN{FS=OFS='${k}'}{print value OFS $0}' | sed "s/^${k}/${k},/" | sed "1 s/${k}/ksize/" | sed "s/similarity/containment/" | sed "s/.prot.sig.gzip//" | sed "s/.cds.name_change.faa//" > ${wd}/compare.prot.${k}.csv

done


