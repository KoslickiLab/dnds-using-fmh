#!/usr/bin/env python
import argparse
import subprocess
import os
import sys
from pathlib import Path
from loguru import logger

"""SKETCH FILES"""

def sketch_dna(fasta, klist, scaled, out_sigfile):
    """sketch a DNA fasta file"""
    """fasta: Fasta input file"""
    """klist: list of k-mer sizes"""
    """scaled: scaled factor"""
    """outfile: signature name"""
    cmd = f"sourmash sketch dna -f -p k={',k='.join(klist)},scaled={scaled} {fasta} -o {out_sigfile} "
    try:
        logger.info(f"Sketching DNA fasta file: {fasta}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully sketched")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while sketching {fasta}: {e}")

def sketch_protein(fasta, klist, scaled, out_sigfile):
    """sketch a protein fasta file"""
    """fasta: Fasta input file"""
    """klist: list of k-mer sizes"""
    """scaled: scaled factor"""
    """outfile: signature name"""
    cmd = f"sourmash sketch protein -f -p k={',k='.join(klist)},scaled={scaled} {fasta} -o {out_sigfile}"
    try:
        logger.info(f"Sketching DNA fasta file: {fasta}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully sketched")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while sketching {fasta}: {e}")

"""CONCATENATE SIGNATURES"""

def cat_signatures(signatures_dir, out_cat_sigfile):
    """Concatenate signature files. This runs when a user has a list of genome fasta files"""
    """signature_list: list of signatures to concatenate"""
    """outfile: signature name"""
    cmd = f"mkdir signatures_dir/../signature_cat"
    cmd = f"sourmash signature cat {signatures_dir}/* -o ${signatures_dir}/../signature_cat/{out_cat_sigfile}"
    try:
        logger.info(f"Making a new directory in ${signatures_dir}../ called signature_cat")
        subprocess.run(cmd, shell=True, check=True)
        logger.info(f"Concatenating individual genome signatures in: {signatures_dir}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully concatenated")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while concatenating for signatures in {signatures_dir}: {e}")

"""COMPARE SIGNATURES VIA CONTAINMENT INDEX"""

def create_compare_directory(dir,content):
    """directory creation to store sourmash compare records"""
    """dir: working directory"""
    """content: identify whether it is dna or protein"""
    cmd = f'mkdir {dir}/compare_{content}'
    try:
        logger.info(f"Creating compare_{content} in {dir}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successful make directory")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while creating directory in {dir}")

def compare_dna_signatures(ref_dna, query_dna, ksize, dir, out_compare_csv):
    """compare dna signature file of ref and query."""
    """ref_dna: reference dna signature"""
    """query_dna: query dna signature"""
    """outfile: signature name"""
    cmd = f"sourmash compare {ref_dna} {query_dna} --containment --dna --ksize {ksize} --csv {dir}/{out_compare_csv}"
    try:
        logger.info(f"Comparing and obtaining containment index between {ref_dna} and {query_dna}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully compared")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while comparing {ref_dna} and {query_dna}: {e}")

def compare_protein_signatures(ref_protein, query_protein, ksize, dir, out_compare_csv):
    """compare dna signature file of ref and query."""
    """ref_dna: reference dna signature"""
    """query_dna: query dna signature"""
    """outfile: signature name"""
    cmd = f"sourmash compare {ref_protein} {query_protein} --containment --protein --ksize {ksize} --csv {dir}/{out_compare_csv}"
    try:
        logger.info(f"Comparing and obtaining containment index between {ref_protein} and {query_protein}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully compared")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while comparing {ref_protein} and {query_protein}: {e}")



