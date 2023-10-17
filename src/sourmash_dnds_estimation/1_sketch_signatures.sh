#!/bin/bash
set -eoux pipefail

wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/ecoli_10_strains_pairwise_genome_sketches/test_small'
format='cds.fna'

cd ${wd}
mkdir signatures

filelist=*.name_change.${format}
echo ${filelist}
for file in $filelist; do
    #echo $file
    SUBSTRING1=$(echo $file | cut -d'.' -f 1,2)
    #other parameters
    scaled=1 #scale factor for reference

    #fasta files
    ref=${file}
    #signature output filenames
    ref_output=signatures/${SUBSTRING1}.dna.sig.gzip

    #fasta files
    sample=${SUBSTRING1}.name_change.cds.faa
    #signature output filenames
    samp_output=signatures/${SUBSTRING1}.prot.sig.gzip

    #sourmash sketch dna -p k=15,k=18,k=21,k=24,k=27,k=30,k=33,k=36,k=39,k=42,k=45,k=60,k=75,k=90,scaled=$scaled $ref -o $ref_output

    nohup sourmash sketch dna -p k=15,k=21,k=30,k=45,k=60,scaled=$scaled $ref -o $ref_output 2>&1 &
    #sourmash sketch protein -p k=5,k=6,k=7,k=8,k=9,k=10,k=11,k=12,k=13,k=14,k=15,k=20,k=25,k=30,scaled=$scaled $samples -o $samp_output
    nohup sourmash sketch protein -p k=5,k=7,k=10,k=15,k=20,scaled=$scaled $sample -o $samp_output 2>&1 &

 
done

