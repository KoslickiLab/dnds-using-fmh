#!/usr/bin/env python

"""SOURMASH BRANCHWATER SCRIPTS""" 

import subprocess
from loguru import logger
import time

def run_manysketch(fasta_file_csv,ksize,scaled,cores, working_dir,molecule=None):
    """sketch multiple dna or protein signature file of ref and query.
    fasta_file_csv: csv file that contains three columns (name, genome_filename, protein_filename) as discussed in sourmash branchwater
    kaa: k-mer size at an amino acid level
    scaled: scaled factor
    molecule: identify the DNA or protein, ksizes depend on molecule
    cores:
    """
    if molecule not in ("dna", "protein", None):
        raise ValueError(f"molecule must be 'dna' or 'protein', got {molecule!r}")
    if molecule=="dna":
        cmd = f'sourmash scripts manysketch {fasta_file_csv} -p k={ksize},scaled={scaled} -c {cores} -o {working_dir}/dna.zip'
    else:
        cmd = f'sourmash scripts manysketch {fasta_file_csv} -p DNA,k={ksize*3},scaled={scaled} -p protein,k={ksize},scaled={scaled} -c {cores} -o {working_dir}/data.zip'
    try:
        logger.info(f"Sketching data fasta file: {fasta_file_csv}")
        start_time = time.time()
        logger.info(f"{cmd}")
        subprocess.run(cmd, shell=True, check=True)
        time_logged = time.time()-start_time
        logger.success(f"Successfully sketched data in {time_logged} seconds")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while sketching {fasta_file_csv}: {e}")
        raise
    


def run_pairwise(zipfile, k, scaled, out_csv, cores, molecule, threshold):
    """compare dna or protein signature file of ref and query
    ref_zipfile: reference dna or protein signature zip file that was produced in manysketch
    query_zipfile: reference dna or protein signature zip file that was produced in manysketch
    kaa: k-mer size at an amino acid level
    scaled: scaled factor
    molecule: identify the DNA or protein, ksizes depend on molecule
    cores:
    working_dir: working directory where to output results"""
    cmd = f"sourmash scripts pairwise {zipfile} -k {k} -s {scaled} -m {molecule} -t {threshold} -o {out_csv} --cores {cores}"
    try:
        logger.info(f"Obtaining pairwise containment index for {zipfile}")
        logger.info(f"{cmd}")
        start_time = time.time()
        subprocess.run(cmd, shell=True, check=True)
        time_logged = time.time()-start_time
        logger.success(f"Successfully pairwise compared in {time_logged} seconds")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while comparing {zipfile}: {e}")
        raise



