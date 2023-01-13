"""Module containing code for predicting open reading frames"""

from more_itertools import sliced
import subprocess
import pandas as pd

def sketch(inputfile1, inputfile2, kmer_size=7,scaledfile1=100, scaledfile2=1):
    """Function that skatches signatures for refernce genome and sample"""
    cmd1=f"sourmash sketch protein -p k={kmer_size}, scaled={scaledfile1} {inputfile1} -o ref-genome.sig"
    subprocess.run(cmd1, stdout=subprocess.PIPE, shell=True)

    cmd2=f"sourmash sketch protein -p k={kmer_size}, scaled={scaledfile2} {inputfile2} -o samples.sig.zip"
    subprocess.run(cmd2, stdout=subprocess.PIPE, shell=True)

def prefetch(sig1, sig2, kmer_size=7, Tbp=1, OUTPUT_FILENAME='prefetch_res,csv'):
    """Function that create MANIFEST CSV file"""
    cmd=f"sourmash prefetch {sig1} {sig2} --protein -o {OUTPUT_FILENAME} --threshold-bp {Tbp} -k {kmer_size}"
    subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)

def search(MANIFEST_CSV_FILE, sig1, sig2):
    """Function that reports containment index of signature files"""
    with open(MANIFEST_CSV_FILE, 'r', encoding='utf-8') as file:
     lines = file.readlines()
     for line in lines:
         if line[0] != "#" and line[0]!='':
             md5_temp = line.split(',')[1]
             if md5_temp!='md5':
                 cmd = f"sourmash search --md5 {md5_temp} {sig2} {sig1} --containment -o {md5_temp}.csv"
                 subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)

def predict_ORF(md5_csv):
    """Function that predicts ORF"""
    data = pd.read_csv({md5_csv}, sep=",", header=None, engine="python", on_bad_lines="skip", quoting=3)
    



