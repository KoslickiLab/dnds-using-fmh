# FRACdNdS

FRACdNdS is a scalable dN/dS estimator for genomes that leverages the sketching approach FracMinHash.

FRACdNdS estimates dN/dS at a genomic level by deriving the containment indexes for between two genomes at a nucleotide and protein level. In practice, the dN/dS ratio infers selection pressues of protein-coding genes by estimating the rates of nonsynonymous substitutions to synonymous ones. We developed a dN/dS estimator to apply on pairwise analyses of longer sequences such as genomes. Currently, FRACdNdS is applied on a FASTA files and has two modes of execution: (1) "sngl" for using a single FASTA with multiple genome entries and (2) "bwpair" for using multiple individual FASTA files, where each FASTA file contains the gene sequences of one genome. 

For detailed usage and installation instructions please visit: coming soon.

Citing FRACdNdS

Please cite: coming soon.

# Requirements

FRACdNdS >= v is written in python3 and requires the following dependencies:

- X
- Y
- Z

# Installation

## Conda environment

under construction

```
# Create the conda ennvironment
conda create -n fracdnds

# Activate the conda environment
conda activate fracdnds

# Download and install dependencies
conda install -c bioconda -c conda-forge XYZ

# Clone the github directory
git clone XYZ
cd FRACdNdS
python setup.py install
```

# Quick start

demo coming soon

To run, simply invoke the `nohup bash run_fmh_omega.sh > log 2>&1 &` script. Other scripts are supporting scripts.

If you need memory information, invoke `nohup /usr/bin/time -v bash run_fmh_omega.sh > log 2>&1 &`

# Usage 

## Overview

## Workflow

