# dn-ds-using-fracminhash
Calculate dn/ds ratio, which is the ratio of nonsynonymous mutations to synonymous mutations, using FracMinHash sketches.

All mathematical calculations are performed here: https://www.dropbox.com/s/gbw4unbs4ewz7e9/DnDs.nb?dl=0 by David.

## Installation

```
conda create -y --name dnds
conda install -y --name dnds -c conda-forge -c bioconda --file requirements.txt
conda activate dnds
```
