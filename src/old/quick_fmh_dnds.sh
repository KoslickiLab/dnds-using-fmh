#!/bin/bash
set -eoux pipefail

wd='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/fmh_dnds_not_using_sourmash_translate/negative_selection'

dna=${wd}/'sequences.fna'
protein=${wd}/'sequences.faa'

echo Sketching of FASTA Files has begun...

#sourmash sketch singleton
dna_sig=${wd}/signatures/sequences.dna.sig.gzip
query_dna_sig=${wd}/signatures/sequences.dna.sig.gzip

protein_sig=${wd}/signatures/sequences.prot.sig.gzip
query_protein_sig=${wd}/signatures/sequences.prot.sig.gzip
scaled=1

echo Creating signatures folder.
mkdir ${wd}/signatures

echo Running sourmash sketch dna using $dna
time sourmash sketch dna -p k=15,k=21,scaled=$scaled $dna --singleton -o $dna_sig > sketchdna.log 2>&1 &
echo Produced the following signature file: $dna_sig

echo Running sourmash sketch protein using $protein
time sourmash sketch protein -p k=5,k=7,scaled=$scaled $protein --singleton -o $protein_sig >sketchprot.log 2>&1 & 
echo Produced the following signature file: $protein_sig

#sourmash compare
echo Creating compare_dna folder.
mkdir ${wd}/compare_dna
echo Creating compare_protein folder.
mkdir ${wd}/compare_protein
for ksize in 5 7; do

    declare -i nt_k=ksize*3
    echo ksize
    echo nt_k
    echo Running sourmash compare on $dna_sig
    echo time sourmash compare $dna_sig $query_dna_sig --containment --dna --ksize $nt_k --o ${wd}/compare_dna/compare.dna.${nt_k}.mat --csv ${wd}/compare_dna/compare.dna.${nt_k}.csv > ${wd}/compare_dna/compare.dna.${nt_k}.log 2>&1 & 
    time sourmash compare $dna_sig $query_dna_sig --containment --dna --ksize $nt_k --o ${wd}/compare_dna/compare.dna.${nt_k}.mat --csv ${wd}/compare_dna/compare.dna.${nt_k}.csv > ${wd}/compare_dna/compare.dna.${nt_k}.log 2>&1 & 
    echo Produced containment files from $dna_sig

    echo Running sourmash compare on $protein_sig
    time sourmash compare $protein_sig $query_protein_sig --containment --protein --ksize $ksize --o ${wd}/compare_protein/compare.prot.${ksize}.mat --csv ${wd}/compare_protein/compare.prot.${ksize}.csv > ${wd}/compare_protein/compare.prot.${ksize}.log 2>&1 & 
    echo Produced containment files from $protein_sig

done

#estimate fmh dnds
for k in 5 7
do
declare -i nt_k=k*3

nt_containment=${wd}/compare_dna/compare.dna.${nt_k}.csv
protein_containment=${wd}/compare_protein/compare.prot.${k}.csv

out=dnds_constant_${k}.csv

echo Estimate FMH dNdS k=${k} from compare_dna and compare_protein found in $wd 
python3 dnds_CI.py --nt ${nt_containment} --protein ${protein_containment} --k ${k} --o ${out} --wd ${wd}/

done

#echo ${wd}dnds_constant_*
echo FMH dNdS estimations done
