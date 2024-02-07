#!/bin/bash
set -eoux pipefail

wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/gKaKs_paper/data'
#file='positive_selection_queries_10002_0.001'
#format='fna'

mkdir ../signatures
cd ${wd}
filelist=*
echo ${filelist}

for file in $filelist; do
    SUBSTRING1=$(echo $file | cut -d'.' -f 1)

    #sourmash sketch parameters
    scaled=1 #scale factor for reference

    #fasta files
    ref=${file}
    #signature output filename for reference
    ref_output=${wd}/../signatures/${SUBSTRING1}.dna.sig.gzip
    #signature output filename for query
    samp_output=${wd}/../signatures/${SUBSTRING1}.prot.sig.gzip

    #sourmash sketch
    nohup time sourmash sketch dna -p k=15,k=21,k=30,k=36,scaled=$scaled $ref -o $ref_output 2>&1 &
    nohup time sourmash sketch protein -p k=5,k=7,k=10,k=12,scaled=$scaled $ref -o $samp_output 2>&1 &

    #sourmash sketch singleton
    #nohup time sourmash sketch dna -p k=15,k=21,k=30,k=45,k=60,scaled=$scaled $ref --singleton -o $ref_output 2>&1 &
    #nohup time sourmash sketch translate -p k=5,k=7,k=10,k=15,k=20,scaled=$scaled $ref --singleton -o $samp_output 2>&1 &
 
done
