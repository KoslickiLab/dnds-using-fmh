#!/bin/sh
"""
GenomegaMap estimates a dN/dS constant among ortholog genes for within species strains. 
Due to this method, sequences where we want to esitmate dN/dS individually to a reference sequence has to be done from individual fasta files in produce dN/dS estimates from GenomegaMap. 
Further, XML files need to be created for each of these fasta files as well.
This bash script is used to create these XML files in order to run GenomegaMap.
"""

#mkdir XML_files_for_sequences_of_10002nt

DIRECTORY=/data/jzr5814/software/genomegaMap/example_10002nt/XML_files_for_sequences_of_10002nt/
END=99
for ((i=0;i<=END;i++)); do
echo $i
    cat /data/jzr5814/software/genomegaMap/examples/porB3.genomegaMapConstantML.a.xml | sed 's/porB3.carriage.noindels.txt/fasta_files_for_sequences_of_10002nt\/query_gene_0.001_'$i'.fna/' | sed 's/libgenomegaMap.so/..\/libgenomegaMap.so/' > $DIRECTORY/10002nt_query_gene_0.001_$i.genomegaMapConstantML.a.xml
    done
