#!/bin/bash
set -eoux pipefail

wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/marine_bacteria/'
format='fna'


cd ${wd}
mkdir signatures

filelist=*.${format}
echo ${filelist}
for file in $filelist; do
    #echo $file
    #SUBSTRING1=$(echo $file | cut -d'.' -f 1,2)
    SUBSTRING1=$(echo $file | cut -d'.' -f 1)
    #other parameters
    scaled=1 #scale factor for reference

    #fasta files
    ref=${file}
    #signature output filenames
    ref_output=signatures/${SUBSTRING1}.dna.sig.gzip

    #fasta files
    #sample=${SUBSTRING1}.cds.name_change.faa
    #signature output filenames
    samp_output=signatures/${SUBSTRING1}.prot.sig.gzip

    #sourmash sketch
    #nohup time sourmash sketch dna -p k=15,k=21,k=30,k=45,k=60,scaled=$scaled $ref -o $ref_output 2>&1 &
    #nohup time sourmash sketch translate -p k=5,k=7,k=10,k=15,k=20,scaled=$scaled $ref -o $samp_output 2>&1 &

    #sourmash sketch singleton
    nohup time sourmash sketch dna -p k=15,k=21,k=30,k=45,k=60,scaled=$scaled $ref --singleton -o $ref_output 2>&1 &
    nohup time sourmash sketch translate -p k=5,k=7,k=10,k=15,k=20,scaled=$scaled $ref --singleton -o $samp_output 2>&1 &
 
done
