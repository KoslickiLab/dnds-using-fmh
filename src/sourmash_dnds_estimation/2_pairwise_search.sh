
#!/bin/bash
time -v
set -eoux pipefail

start=`date +%s`

wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_genome_sketches'
format='sig.gzip'

cd ${wd}
mkdir compare_protein
mkdir compare_dna

output_compare_protein=${wd}/compare_protein
output_compare_dna=${wd}/compare_dna

cd ${wd}/signatures
filelist=*.${format}

for file1 in $filelist; do
    for file2 in $filelist; do
        SUBSTRING1=$(echo $file1 | cut -d'.' -f 1,2)
        SUBSTRING2=$(echo $file2 | cut -d'.' -f 1,2)
        echo $SUBSTRING1
        echo $SUBSTRING2
        if [ ${SUBSTRING1} != ${SUBSTRING2} ]; then
            for ksize in 5 7 10 15 20; do
                nt_k=$(( 3*ksize ))
                echo $ksize
                #if ! test -f compare_dna/compare.dna.${nt_k}.${SUBSTRING2}_${SUBSTRING1}.csv; then
                nohup time sourmash search $SUBSTRING1.dna.sig.gzip $SUBSTRING2.dna.sig.gzip --containment --dna --ksize $nt_k --t 0 --o ${output_compare_dna}/compare.dna.${nt_k}.${SUBSTRING1}_${SUBSTRING2}.csv > ${output_compare_dna}/compare.dna.${nt_k}.${SUBSTRING1}_${SUBSTRING2}.txt 2>&1 & #assign to diff log files
                #fi
                #if ! test -f compare_protein/compare.prot.${ksize}.${SUBSTRING2}_${SUBSTRING1}.csv; then
                nohup time sourmash search $SUBSTRING1.prot.sig.gzip $SUBSTRING2.prot.sig.gzip --containment --protein --ksize $ksize --t 0 --o ${output_compare_protein}/compare.prot.${ksize}.${SUBSTRING1}_${SUBSTRING2}.csv > ${output_compare_protein}/compare.prot.${ksize}.${SUBSTRING1}_${SUBSTRING2}.txt 2>&1 & #assign to diff log files
                #fi
            done
        fi
    done
done


time -v