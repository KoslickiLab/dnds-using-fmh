

# create axt file from a fasta file where the first sequence is the reference
awk '$0 ~ ">" {c=substr($0,2,length($0))} NR == 2 {ref=$0} NR % 4 == 0 {print c"\n"ref"\n"$0"\n"}' sequences.fna > kaks_sequences.axt

# run kaks_calculator on axt file
KaKs_Calculator -i /data/jzr5814/kaks_calc_tool_analysis/real_data_0.001_10002/kaks_sequences.axt -o /data/jzr5814/kaks_calc_tool_analysis/real_data_0.001_10002/kaks_sequences.axt.kaks -m NG -m YN -m MYN -m gNG -m gYN -m gMYN




