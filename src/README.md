# dn-ds-using-containment

Calculate dn/ds ratio using the containment index. This program takes as input protein sequences.

#The command to reproduce tests is shown below:

# Reproduce tests
`python3 dnds_using_containment_index.py --fasta1 uniprotkb50.fasta --fasta2 50_sequences.fna --predict orf --k 7,14 --wd results/`

To run in the background use: `nohup python3 dnds_using_containment_index.py --fasta1 uniprotkb50.fasta --fasta2 50_sequences.fna --predict orf --k 7,14 --wd results/ > my_log.txt 2>&1 &`

For additional help and parameters run `python3 dnds_using_containment_index.py --h`
