# dn-ds-using-containment

Calculate dn/ds ratio using the containment index. This program takes as input protein sequences.

#The command to reproduce tests is shown below:

# Reproduce tests
`python3 dnds_using_containment_index.py --ref_fasta ref.fna --query_fasta query.fna --k 3,4,5,6 --wd WORKING_DIRECTORY_FOLDER --analyze dna`

To run in the background use: `nohup python3 dnds_using_containment_index.py --ref_fasta ref.fna --query_fasta query.fna --k 3,4,5,6 --wd WORKING_DIRECTORY_FOLDER --analyze dna > nt_log.txt 2>&1 &`

For additional help and parameters run `python3 dnds_using_containment_index.py --h`
