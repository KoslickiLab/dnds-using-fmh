"""Module containing code for predicting open reading frames"""

from more_itertools import sliced
import subprocess
import pandas as pd

def get_kmer_argument(kmer_list):
    klst = kmer_list.split(',')
    kmers=[]
    for k in klst:
        kmers.append("k="+str(k))
    if len(kmers) > 1:
        kmers_arg = ",".join(kmers)
    else:
        kmers_arg = kmers[0]
    return(kmers_arg)

def sketch(inputfile1, inputfile2, kmer_size=7,scaledfile1=100, scaledfile2=1, outputfile1="ref-genome.sig", outputfile2="samples.sig.zip"):
    """Function that skatches signatures for refernce genome and sample"""
    kmer_sizes = get_kmer_argument(kmer_size)

    cmd1=f"sourmash sketch protein -p {kmer_sizes},scaled={scaledfile1} {inputfile1} -o {outputfile1}"
    print(cmd1)
    subprocess.run(cmd1, stdout=subprocess.PIPE, shell=True)

    cmd2=f"sourmash sketch protein -p {kmer_sizes},scaled={scaledfile2} {inputfile2} --singleton -o {outputfile2}"
    subprocess.run(cmd2, stdout=subprocess.PIPE, shell=True)

def prefetch(sig1, sig2, kmer_size=7, Tbp=1, OUTPUT_FILENAME='prefetch_res.csv'):
    """Function that create MANIFEST CSV file"""
    kmer_sizes = get_kmer_argument(kmer_size).replace("k=","").split(",")
        
    if len(kmer_sizes) > 1:
        with open("/home/grads/jzr5814/data/prefetch.txt", 'w', encoding="utf-8") as output:
            for kmer in kmer_sizes:
                cmd1=f"sourmash prefetch {sig1} {sig2} --protein -o /home/grads/jzr5814/data/prefetch_res_{kmer}.csv --threshold-bp {Tbp} -k {kmer}" 
                output.write(cmd1+"\n")
            output.close()
        cmd2 = f"parallel -d '\n' < /home/grads/jzr5814/data/prefetch.txt"
        subprocess.run(cmd2, stdout=subprocess.PIPE, shell=True)
    else:
        cmd=f"sourmash prefetch {sig1} {sig2} --protein -o /home/grads/jzr5814/data/prefetch_res.csv --threshold-bp {Tbp} -k {kmer_size}"
        subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)

def search(MANIFEST_CSV_FILE, sig1, sig2, output_directory):
    """Function that reports containment index of signature files"""

    with open(MANIFEST_CSV_FILE, 'r', encoding='utf-8') as file:
     lines = file.readlines()
     for line in lines:
         if line[0] != "#" and line[0]!='':
             md5_temp = line.split(',')[1]
             if md5_temp!='md5':
                 cmd = f"sourmash search --md5 {md5_temp} {sig2} {sig1} --containment -o {output_directory}/{md5_temp}.csv"
                 subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)

def predict_ORF(md5_csv):
    """Function that predicts ORF"""
    data = pd.read_csv({md5_csv}, sep=",", header=None, engine="python", on_bad_lines="skip", quoting=3)
    



