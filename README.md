
![example workflow](https://github.com/Dreycey/DebruijnExtend/actions/workflows/github_actions.yml/badge.svg)

![Debruijn Extend](figures/debruijnextend_logo.png)

# DebruijnExtend
### This repository contains the source code for the DebruijnExtend tool. Overall, this tool predicts the secondary structure of a protein(s) given an input fasta file with protein(s).

DebruijnExtend is a machine learning algorithm that constructs a probability-based de Bruijn graph using kmers of protein secondary structure fragments associated with a given primary sequence. Once created, edge contraction is used to find paths of highest probability to determine the overall predicted secondary structure.

# Installation
USING CONDA:
```
conda env create -f environment.yml
```

# How to use this tool
The most detailed information on usage is featured on the [DebruijnExtend Wiki](https://github.com/Dreycey/DebruijnExtend/wiki)

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

## Basic Usage (Building Data Locally)

### cluster-heuristic KNN (fast) (not using default file)
EXAMPLE:                                                                        
```
python3 DebruijnExtend.py -i examples/gfp.fasta \
                          -o result.ss3
                          -ht 6 \
                          -c 4
```

**note**: This first has to create the hash tables and KNN clusters, so it may take a couple minutes to create the needed information.

## Basic Usage (Using downloaded Default Data)

### Download Defualt files
Use a tool, `gdown`, to retrieve the files. `gdown` is included in the conda file.
```
gdown --folder --id 1VJKvfBaNFQj0zXCwKDTKcBzkgRT26bz-
```

### cluster-heuristic KNN (using default files downloaded using gdown)
EXAMPLE:                                                                        
```
python3 DebruijnExtend.py -i examples/gfp.fasta \
                          -o result.ss3
                          -ht DefaultData/prothashtable_10.p \
                          -c DefaultData/cluster_file.pickle
```
**note**: Before using this method, make sure to get the premade hash table and clusters using the gdown command ('Download Defualt files' section above)
