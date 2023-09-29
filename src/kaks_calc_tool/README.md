# if fasta names are compliated in which they have spaces and other chars use the following command to add quotations (i.e. you may need to use this for seq.parse fasta dictionary)
cat sequences_with_no_line_breaks.fna | sed 's/>\(.*\)/>"\1"/g' | grep '>'

# fasta sequences shouldn't have new lines in order for the program to give you correct results
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' sequence_temp.fna > sequence.fna

rm sequence_temp

# create axt file from concatenated alignment fasta files
# sed '0~4G' inserts a new line at every 4 line
# sed '3~5d' removes every 5th line starting at the third line
# sed 's/-$/--/g' adds an extra gap in case of any sequences that cannot be divided by 3
cat clustalw_DNA_all.aln | sed '0~4G'| sed '3~5d' > kaks_sequences.axt

# removes new line and parses a fasta file into an axt file
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' | sed '3~4d' | sed '0~3 a\\'

# run kaks_calculator on axt file
KaKs_Calculator -i kaks_sequences.axt -o kaks_sequences.axt.kaks -m GNG -m GY -m LPB -m LWL -m NG -m YN




