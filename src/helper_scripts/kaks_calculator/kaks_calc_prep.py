from Bio.SeqUtils.CheckSum import seguid
from Bio import SeqIO

'''
This script is to create a AXT file meant for pairwise dn/ds estimates using the kaks_calculator.

If you have a fasta file with multiple entries, then you want to use this script to create an axt file for pairwise dn/ds.

This script was mean't to run on a fasta file generated with a random sequence of length X. The length has to be divisible by 3 or the kaks_calculator will return an error.

If you have sequences of differing lengths, then using an aligner such as clustalw can help but another script would have to be implemented to create a AXT file for kaks_calculator input. 
'''

input_file = 'negative_selection_queries_10002_0.01.fna'
output = open('negative_selection_queries_10002_0.01.fna.axt','w')

seguid_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"), lambda rec: seguid(rec.seq))

for key in seguid_dict:
    ref = seguid_dict[key].description
    ref_seq = str(seguid_dict[key].seq)
    for key in seguid_dict:
        output.write(ref)
        output.write('\n')
        output.write(ref_seq)
        output.write('\n')
        output.write(str(seguid_dict[key].seq))
        output.write('\n')
        output.write('\n')

output.close()