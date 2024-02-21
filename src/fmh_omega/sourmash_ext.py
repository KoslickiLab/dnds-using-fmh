#!/usr/bin/env python
import subprocess
from loguru import logger

"""SKETCH FILES"""
def sketch_genome_dna(fasta, klist, scaled, out_sigfile, multiple=None):
    """sketch a DNA fasta file. Used when a fasta include genomic information of one species
    fasta: Fasta input file
    klist: list of k-mer sizes
    scaled: scaled factor
    outfile: signature name"""
    if multiple: #Used when a fastn include genomic information of one species
        cmd = f"sourmash sketch dna -f -p k={klist},scaled={scaled} {fasta} -o {out_sigfile} "
    else: #Used when a fastn file include multiple species entries
        cmd = f"sourmash sketch dna -f -p k={klist},scaled={scaled} --singleton {fasta} -o {out_sigfile} "
    try:
        logger.info(f"Sketching DNA fasta file: {fasta}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully sketched DNA")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while sketching {fasta}: {e}")

def sketch_genome_protein(fasta, klist, scaled, out_sigfile, multiple=None):
    """sketch a protein fasta file. 
    fasta: Fasta input file
    klist: list of k-mer sizes
    scaled: scaled factor
    outfile: signature name"""
    if multiple: #Used when a fastn include genomic information of one species
        cmd = f"sourmash sketch protein -f -p k={klist},scaled={scaled} {fasta} -o {out_sigfile} "    
    else: #Used when a fastn file include multiple species entries.
        cmd = f"sourmash sketch protein -f -p k={klist},scaled={scaled} --singleton {fasta} -o {out_sigfile}"
    try:
        logger.info(f"Sketching DNA fasta file: {fasta}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully sketched protein")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while sketching {fasta}: {e}")

"""COMPARE SIGNATURES VIA CONTAINMENT INDEX"""
def compare_signatures(ref, query, ksize, molecule, working_dir):
    """compare dna or protein signature file of ref and query."""
    """ref: reference dna or protein signature"""
    """query: query dna or protein signature"""
    """outfile: signature name"""
    cmd = f"sourmash compare {ref} {query} --containment --{molecule} --ksize {ksize} --csv {working_dir}/compare_{molecule}/compare.{molecule}.{ksize}.csv"
    try:
        logger.info(f"Comparing and obtaining containment index between {ref} and {query}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully compared")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while comparing {ref} and {query}: {e}")

"""SOURMASH BRANCHWATER SCRIPTS"""
def run_manysketch(fasta_file_csv,klist,scaled,cores,molecule, working_dir):
    """sketch multiple dna or protein signature file of ref and query.
    fasta_file_csv: csv file that contains three columns (name, genome_filename, protein_filename) as discussed in sourmash branchwater
    klist: list of k-mer sizes
    scaled: scaled factor
    molecule: identify the list of ksizes, ksizes depend on molecule
    cores:
    """
    cmd = f'sourmash scripts manysketch {fasta_file_csv} -p {molecule},k={klist},scaled={scaled} -c {cores} -o {working_dir}/{molecule}.zip'
    try:
        logger.info(f"Sketching {molecule} fasta file: {fasta_file_csv}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully sketched {molecule}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while sketching {fasta_file_csv}: {e}")

def run_multisearch(ref_zipfile,query_zipfile, ksize, scaled, out_csv, cores, molecule):
    """compare dna or protein signature file of ref and query
    ref_zipfile: reference dna or protein signature zip file that was produced in manysketch
    query_zipfile: reference dna or protein signature zip file that was produced in manysketch
    klist: list of k-mer sizes
    scaled: scaled factor
    molecule: identify the list of ksizes, ksizes depend on molecule
    cores:
    working_dir: working directory where to output results"""
    cmd = f"sourmash scripts multisearch {ref_zipfile} {query_zipfile} -k {ksize} -s {scaled} -m {molecule} -t 0 -o {out_csv} --cores {cores}"
    try:
        logger.info(f"Comparing and obtaining containment index between {ref_zipfile} and {query_zipfile}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully compared")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while comparing {ref_zipfile} and {query_zipfile}: {e}")

def run_pairwise(zipfile, ksize, scaled, out_csv, cores, molecule):
    """compare dna or protein signature file of ref and query
    ref_zipfile: reference dna or protein signature zip file that was produced in manysketch
    query_zipfile: reference dna or protein signature zip file that was produced in manysketch
    klist: list of k-mer sizes
    scaled: scaled factor
    molecule: identify the list of ksizes, ksizes depend on molecule
    cores:
    working_dir: working directory where to output results"""
    cmd = f"sourmash scripts pairwise {zipfile} -k {ksize} -s {scaled} -m {molecule} -t 0 -o {out_csv} --cores {cores}"
    try:
        logger.info(f"Obtaining pairwise containment index for {zipfile}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully compared")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while comparing {zipfile}: {e}")