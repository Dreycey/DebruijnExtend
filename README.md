
![example workflow](https://github.com/Dreycey/DebruijnExtend/actions/workflows/github_actions.yml/badge.svg)

![Debruijn Extend](figures/debruijnextend_logo.png)

# DebruijnExtend
### This repository contains the source code for the DebruijnExtend tool. Overall, this tool predicts the secondary structure of a protein given an input fasta file with protein(s).

DebruijnExtend is a machine learning algorithm that constructs a probability-based de Bruijn graph using kmers of protein secondary structure fragments associated with a given primary sequence. Once created, edge contraction is used to find paths of highest probability to determine the overall predicted secondary structure.

# Installation
USING CONDA:
```
conda env create -f environment.yml
```

# Download Defualt files
Use a tool, `gdown`, to retrieve the files. `gdown` is included in the conda file.
```
gdown --folder --id 1VJKvfBaNFQj0zXCwKDTKcBzkgRT26bz-
```

# How to use this tool

## General Usage
* General Usage
```
DebruijnExtend.py [-h] [-v] -i INPUT (protein fasta) \
                            -o OUTPUT_FILE (ss3 file) \
                            -ht HASH_TABLE \
                            [-c USE_CLUSTERS] \
                            [-d CSV_DATA_PATH] \
                            [-t THREADS]
```

The input fasta file may be a single protein or a multifasta. The ordering of the output file will have the same format with the predicted 3-state secondary structure appearing below the name of the protein.

## Using DefaultData (pulled using the above)
Before using the below methods, make sure to get the premade hash table and clusters using the gdown command ('Download Defualt files' section above)

### using default KNN (non-cluster method; time expensive)

EXAMPLE (kmer size=10,CT=6):                                                                        
```
python3 DebruijnExtend.py -i examples/gfp.fasta \
                          -o result.ss3 \
                          -ht DefaultData/prothashtable_10.p
```

### using default file with cluster-heuristic KNN (fast)
EXAMPLE:                                                                        
```
python3 DebruijnExtend.py -i examples/gfp.fasta \
                          -o result.ss3
                          -ht DefaultData/prothashtable_10.p \
                          -c DefaultData/cluster_file.pickle
```

## Creating new data files (hash table and clusters)

### Creating new hash table (with default CSV)
If the goal is to create a hash table of size 6, an integer can be supplied to the kmer flag `-ht 6`:
```
Python DebruijnExtend.py -i examples/gfp.fasta \
                         -o result.ss3 \
                         -ht 6
```

### Creating new hash table AND clusters (with default CSV)
If the goal is to create a hash table of size 6, an integer can be supplied to the kmer flag `-ht 6`:
```
Python DebruijnExtend.py -i examples/gfp.fasta \
                         -o result.ss3 \
                         -ht 6 \
                         --use_clusters 3
```
This creates a hash table of kmer size 6, if it doesn't exists. This hash table is then used to create the clusters with a hamming distance cut off of 3.

### Creating new hash table (with new CSV)
A new CSV for creating the hash table and clusters can be used, as long as a path is given:
```
Python DebruijnExtend.py -i examples/gfp.fasta \
                         -o result.ss3 \
                         -ht 10 \
                         -c 7 \
                         -d ${PWD}/NEW.csv 
```
 
**note**, any supplied CSV must follow the following format:
```
sequence length,PDB name,Protein Sequence,3 char ss3
50,'1PTQ','HRFKVYNYMSPTFCDHCGSLLWGLVKQGLKCEDCGMNVHHKCREKVANLC','CCEEEECCCCCCECCCCCCECCCCCCCEEEECCCCCEECHHHHCCCCCCC'
```
