#!/usr/bin/env python
import subprocess
from loguru import logger
import time

"""SOURMASH BRANCHWATER SCRIPTS"""
def run_manysketch(fasta_file_csv,kaa,scaled,cores, working_dir,molecule=None):
    """sketch multiple dna or protein signature file of ref and query.
    fasta_file_csv: csv file that contains three columns (name, genome_filename, protein_filename) as discussed in sourmash branchwater
    klist: list of k-mer sizes
    scaled: scaled factor
    molecule: identify the list of ksizes, ksizes depend on molecule
    cores:
    """
    if molecule=="dna":
        cmd = f'sourmash scripts manysketch {fasta_file_csv} -p k={kaa},scaled={scaled} -c {cores} -o {working_dir}/dna.zip'
    else:
        cmd = f'sourmash scripts manysketch {fasta_file_csv} -p DNA,k={kaa*3},scaled={scaled} -p protein,k={kaa},scaled={scaled} -c {cores} -o {working_dir}/data.zip'
    try:
        logger.info(f"Sketching data fasta file: {fasta_file_csv}")
        start_time = time.time()
        logger.info(f"{cmd}")
        subprocess.run(cmd, shell=True, check=True)
        end_time = time.time()
        time_logged = end_time-start_time
        logger.success(f"Successfully sketched data in {time_logged} seconds")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while sketching {fasta_file_csv}: {e}")
    


def run_pairwise(zipfile, k, scaled, out_csv, cores, molecule, threshold):
    """compare dna or protein signature file of ref and query
    ref_zipfile: reference dna or protein signature zip file that was produced in manysketch
    query_zipfile: reference dna or protein signature zip file that was produced in manysketch
    klist: list of k-mer sizes
    scaled: scaled factor
    molecule: identify the list of ksizes, ksizes depend on molecule
    cores:
    working_dir: working directory where to output results"""
    cmd = f"sourmash scripts pairwise {zipfile} -k {k} -s {scaled} -m {molecule} -t {threshold} -o {out_csv} --cores {cores}"
    try:
        logger.info(f"Obtaining pairwise containment index for {zipfile}")
        logger.info(f"{cmd}")
        start_time = time.time()
        subprocess.run(cmd, shell=True, check=True)
        end_time = time.time()
        time_logged = end_time-start_time
        logger.success(f"Successfully pairwise compared in {time_logged} seconds")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while comparing {zipfile}: {e}")



