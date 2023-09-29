#!/bin/bash

# removes new line and parses a fasta file into an axt file
cat /data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise_take_2_with_CDS_transcripts/divisible_by_3.aln | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' | sed '3~4d' | sed '0~3 a\\' > /data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise_take_2_with_CDS_transcripts/kaks.axt

# run kaks_calculator on axt file
nohup KaKs_Calculator -i /data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise_take_2_with_CDS_transcripts/kaks.axt -o /data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise_take_2_with_CDS_transcripts/kaks.axt.kaks -m GNG -m GY -m LPB -m LWL -m NG -m YN > /data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise_take_2_with_CDS_transcripts/kaks.log 2>&1 &


