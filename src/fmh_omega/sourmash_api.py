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
def cat_signatures(working_dir, molecule, out_name):
    """Concatenate signature files. This runs when a user has a list of genome fasta files"""
    """signature_list: list of signatures to concatenate"""
    """molecule: only be dna or protein"""
    cmd = f"sourmash signature cat {working_dir}/signatures/*.{molecule}.sig.gzip -o ${working_dir}/signatures/{out_name}.cat.{molecule}.sig.gzip"
    try:
        logger.info(f"Making a new directory in ${working_dir} called signature_cat")
        subprocess.run(cmd, shell=True, check=True)
        logger.info(f"Concatenating individual genome signatures in: {working_dir}/signatures")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully concatenated")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while concatenating for signatures in {working_dir}/signatures: {e}")

"""COMPARE SIGNATURES VIA CONTAINMENT INDEX"""
def compare_signatures(ref, query, ksize, molecule, working_dir):
    """compare dna or protein signature file of ref and query."""
    """ref: reference dna or protein signature"""
    """query: query dna or protein signature"""
    """outfile: signature name"""
    cmd = f"sourmash compare {ref} {query} --containment --{molecule} --ksize {ksize} --csv {working_dir}/compare_{molecule}/compare.{molecule}.${ksize}.csv"
    try:
        logger.info(f"Comparing and obtaining containment index between {ref} and {query}")
        subprocess.run(cmd, shell=True, check=True)
        logger.success(f"Successfully compared")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred while comparing {ref} and {query}: {e}")

