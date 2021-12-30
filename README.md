
![example workflow](https://github.com/Dreycey/DebruijnExtend/actions/workflows/github_actions.yml/badge.svg)

![Debruijn Extend](figures/debruijnextend_logo.png)

# DebruijnExtend
### This repository contains the source code for the DebruijnExtend tool. Overall, this tool predicts the secondary structure of a protein given an input fasta file with protein(s).

DebruijnExtend is a machine learning algorithm that constructs a probability-based de Bruijn graph using kmers of protein secondary structure fragments associated with a given primary sequence. Once created, edge contraction is used to find paths of highest probability to determine the overall predicted secondary structure.

# Installation
USING CONDA:
```
conda env create -f DebruijnExtend_environment.yaml
```

# How to use this tool
USAGE:
```                                                                           
python DebruijnExtend.py <input fasta> <kmer size> <output file>    
```

EXAMPLE:                                                                        
```
python DebruijnExtend.py examples/gfp.fasta 4 gfp.ss3
```

