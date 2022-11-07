"""Module to shell from python script"""
import subprocess

cmd1=f"sourmash sketch protein -p k=7,k=14,k=21,k=28,k=35,scaled=100 ../data/uniprotkb.fasta -o ../500_sequences_small_kmers/ref-genome.sig"
subprocess.run(cmd1, stdout=subprocess.PIPE, shell=True)

#cmd2=f"sourmash sketch protein -p k=7,k=42,k=49,k=56,k=63,k=77,scaled=100 ../data/uniprotkb.fasta -o ../500_sequences_large_kmers/ref-genome.sig"
#subprocess.run(cmd2, stdout=subprocess.PIPE, shell=True)

#cmd3=f"sourmash sketch protein -p k=7,k=14,k=21,k=28,k=35,scaled=1 ../data/frames_500_sequences.faa --singleton -o ../500_sequences_small_kmers/orfs_500.sig.zip"
#subprocess.run(cmd3, stdout=subprocess.PIPE, shell=True)

#cmd4=f"sourmash sketch protein -p k=7,k=42,k=49,k=56,k=63,k=77,scaled=1 ../data/frames_500_sequences.faa --singleton -o ../500_sequences_large_kmers/orfs_500.sig.zip"
#subprocess.run(cmd4, stdout=subprocess.PIPE, shell=True)



