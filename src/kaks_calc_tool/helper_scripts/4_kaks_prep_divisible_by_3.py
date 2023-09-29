from Bio.SeqUtils.CheckSum import seguid
from Bio import SeqIO
import os

#use alignment fasta files

fasta_file = '/data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise_take_2_with_CDS_transcripts/all.aln'
output_name = '/data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise_take_2_with_CDS_transcripts/divisible_by_3.aln'
output = open(output_name,'w')

#does not for duplicate names checkout readme file
seguid_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"), lambda rec: seguid(rec.seq))


def divisible_by_3(sequence):
    if len(sequence) % 3 != 0:
        if len(sequence+'-') % 3 != 0:
            return(sequence+'--')
        else:
            return(sequence+'-')
    else:
        return(sequence)

for key1 in seguid_dict:
    ref = seguid_dict[key1].description
    ref_seq = str(seguid_dict[key1].seq)
    output.write('>'+ref)
    output.write('\n')
    output.write(divisible_by_3(ref_seq))
    output.write('\n')
output.close()
#cmd=output_name+""" awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' | sed '3~4d' | sed '0~3 a\\' > kaks.axt"""
#os.system(cmd)  # returns the exit code in unix