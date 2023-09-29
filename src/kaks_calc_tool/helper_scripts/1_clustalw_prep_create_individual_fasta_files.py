from Bio.SeqUtils.CheckSum import seguid
from Bio import SeqIO
import os

wd='/data/jzr5814/orthologr_analysis/HIT000324409_pairwise'
format = 'fna'
fasta_file = f'{wd}/sequences_with_no_line_breaks.{format}'

seguid_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"), lambda rec: seguid(rec.seq))

cmd1 = f"mkdir {wd}/individual_fasta_files"
returned_value = os.system(cmd1)  # returns the exit code in unix

for key in seguid_dict:
    id = seguid_dict[key].id
    ref = seguid_dict[key].description
    ref_seq = str(seguid_dict[key].seq)
    species_fasta = f'{wd}/individual_fasta_files/{id}.{format}'
    output = open(species_fasta,'w')
    output.write('>')
    output.write(ref)
    output.write('\n')
    output.write(ref_seq)
    output.write('\n')
    output.close()
    



