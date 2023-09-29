from Bio.SeqUtils.CheckSum import seguid
from Bio import SeqIO


fasta_file = '/data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise/sequences_with_no_line_breaks.fna'
output = open('/data/jzr5814/kaks_calc_tool_analysis/HIT000324409_pairwise/kaks_sequences.axt','w')

seguid_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"), lambda rec: seguid(rec.seq))

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

