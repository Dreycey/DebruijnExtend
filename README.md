# DebruijnExtend
This repository contains the source code for the DebruijnExtend tool. This tool uses known secondary structures for common protein Khmers to predict the most likely secondary structure for a given input primary protein sequence.

# Installation


# How to use this tool
USAGE:
```                                                                           
python DebruijnExtend_v2.py <input fasta> <kmer size> <output file>    
```

EXAMPLE:                                                                        
```
python DebruijnExtend_v2.py examples/gfp.fasta 4 gfp.ss3
```

# Directory Structure

## docs (informational)
This directory contains a document regarding the project, file sizes for the curated csvs, a table outline, and the project presentation.


## presentations (informational)
This directory contains information regarding presentations.

## datasets (DATA)
This dataset directory contains the primary datasets before processing (processing done using the python scripts).

## curated_csvs (DATA)
This directory holds the curated CSV files that may be imported into the SQL table using the "make tables" script int he sql_scripts directory. It contains the primary datasets after processing (processing done using the python scripts).

## sec_struct_dataset (DATA)
This dataset directory contains the primary datasets that pairs protein sequences up with the corresponding secondary structure. These files can be used as a starting point for further developing the algorithms.

## py_scripts (CODE)
This directory contains code for the python files. These were used for dataset preprocessing.

## sql_scripts (CODE)
This directory contains SQL files that can be used to build the databases, as well as predict related GO terms and for predicting secondary structure from primary sequence alone. 

