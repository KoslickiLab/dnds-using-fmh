

## Create dataset.csv file

As sourmash isntructs: create a CSV file with three headers (name, genome_filename, protein_filename)

## Sketch fromfile

```nohup sourmash scripts manysketch datasets.csv -p protein,k=15,k=21,k=27,k=33,scaled=1 -c 4 -o manysketch_dna.zip > manysketch_dna.log 2>&1```

```nohup sourmash scripts manysketch datasets.csv -p protein,k=5,k=7,k=9,k=11,scaled=1 -c 4 -o manysketch_protein.zip > manysketch_protein.log 2>&1```

## Run sourmash scripts multisearch to obtain 

```nohup sourmash scripts multisearch manysketch_dna.zip manysketch_dna.zip -k 15 -s 1 -o results_dna.csv --cores 4 > multisearch_dna.log 2>&1 ```

```nohup sourmash scripts multisearch manysketch_protein.zip manysketch_protein.zip -k 5 -s 1 -o results_protein.csv --cores 4 > multisearch_protein.log 2>&1 ```

## Estimate dNdS