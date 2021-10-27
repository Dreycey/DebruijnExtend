![Debruijn Extend](figures/debruijnextend_logo.png)

# DebruijnExtend
This repository contains the source code for the DebruijnExtend tool. This tool uses known secondary structures for common protein Khmers to predict the most likely secondary structure for a given input primary protein sequence.

# Installation
* Using Conda
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

