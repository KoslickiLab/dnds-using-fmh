

# create axt file from a fasta file where the first sequence is the reference
awk '$0 ~ ">" {c=substr($0,2,length($0))} NR == 2 {ref=$0} NR % 2 == 0 {print c"\n"ref"\n"$0"\n"}' sequences.fna > kaks_sequences.axt

# run kaks_calculator on axt file
KaKs_Calculator -i kaks_sequences.axt -o kaks_sequences.axt.kaks -m NG -m YN -m MYN -m GNG -m GYN -m GMYN -m GY

# KaKs_Caulculator requires that sequences be divisible by three. If you encounter this problem, the following command should fix it.
cat sequences.fna | cut -c 2- | sed 's/^gene/>gene/' | sed 's/^ref/>ref/' > sequences_modified.fna 


