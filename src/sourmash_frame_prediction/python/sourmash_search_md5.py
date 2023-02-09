"""Module to shell from python script"""
import subprocess

with open('../500_sequences/SOURMASH-MANIFEST.csv', 'r', encoding='utf-8') as file:
    lines = file.readlines()
    for line in lines:
        if line[0] != "#" and line[0]!='':
            md5_temp = line.split(',')[1]
            if md5_temp!='md5':
                cmd = f"sourmash search --md5 {md5_temp} ../500_sequences/orfs_500.sig.zip ref-genome.sig --containment -o ../500_sequences/csv_files/{md5_temp}.csv"
                subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
