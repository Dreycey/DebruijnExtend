# prot_hashtables
This directory contains the pickled hash tables used in DebruijnExtend. The hash table is essentially a python dictionary with primary sequences and the associated secondary counts observed in the training CSV input.

Example:
```
hashtable = {"AGCTH" : {"CEEEH" : 2, ..., "EEEEE" : 30},
                ...
             "TTPGH" : {"HHHEE" : 34, ..., "HHEEE" : 2}}
```