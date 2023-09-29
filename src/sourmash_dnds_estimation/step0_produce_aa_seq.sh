#!/bin/bash
set -eoux pipefail

## transseq
# translate nucleotide to amino acid on commandline
transeq /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/HIT000324409_pairwise_take_2_with_CDS_transcripts/pairwise_HIT000324409.fastn /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/HIT000324409_pairwise_take_2_with_CDS_transcripts/pairwise_HIT000324409_temp.fasta

# transeq changes sequence names adding reading frame to the end and if not changed, will return empty fmh dnds
cat /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/HIT000324409_pairwise_take_2_with_CDS_transcripts/pairwise_HIT000324409_temp.fasta | sed 's/_1//' > /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/HIT000324409_pairwise_take_2_with_CDS_transcripts/pairwise_HIT000324409.fasta

rm /data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/HIT000324409_pairwise_take_2_with_CDS_transcripts/pairwise_HIT000324409_temp.fasta
