#! usr/bin/python3
"""
This is the primary script and code for running the DebruijnExtend algorithm.

USAGE
         python DebruijnExtend.py <input fasta> <kmer size> <output file>
EXAMPLE:
         python DebruijnExtend.py examples/gfp.fasta 4 gfp.ss3
         
Updates (v2; 11/28/2020): 
    * Here unelongated sequences are not extended further.
    * There is a max size set for the debruijn extension. This helps
      speed the algorithm dramatically and allows for scaling to big
      protein sequences.

Updates (v3; 12/07/2020): 
    * Here top matches are given random guesses if no edge

Updates (10/27/2021): 
    * refactoring all methods/functions
        * Adding type hinting, imrpoved doc strings
        * getting rid of other "versions"

Updates (10/04/2021): 
    * allowing for multifasta input
"""
# std pkgs
import math
import sys
from typing import Dict, List, Optional, Union
import numpy as np
from pathlib import Path
import os
# import multiprocessing as mp
# import threading
# from concurrent.futures import ThreadPoolExecutor
# import itertools
# non-std pkgs
import random
import argparse
# in-house packages
from src.debruijnextend.KmerCluster import KmerCluster
from src.debruijnextend.DebruijnExtend import DebruijnExtend, HashTableType
from src.debruijnextend.utils import readFasta, get_kmer_size, clean_input_sequences




def parseArgs(argv=None) -> argparse.Namespace:
    """
    This method takes in the arguments from the command and performs
    parsing.
    INPUT: 
        Array of input arguments
    OUTPUT:
        returns a argparse.Namespace object
    """
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-i", "--input", help="the input fasta file", required=True)
    parser.add_argument("-ht", "--hash_table", help="the input hash table (python pickle file)", required=True)
    parser.add_argument("-o", "--output_file", help="the output file prefix", required=True)
    parser.add_argument("-t", "--threads", help="number of threads", required=False)
    parser.add_argument("-c", "--use_clusters", default=False, help="use clustering instead of O(N) all vs all", action='store_true', required=False)
    return parser.parse_args(argv)

# RUN SCRIPT
def main():
    """
    This function controls the flow of the script.
    """
    args = parseArgs(sys.argv[1:])
    # get input ready
    proteins, protein_names = readFasta(args.input)
    hash_table = Path(args.hash_table)
    outputfile = Path(args.output_file)
    # get protein clusters
    if args.use_clusters:
        kmer_clust = KmerCluster.init_struct(hash_table, "cluster_file.pickle")
    # delete output file if exists
    if os.path.exists(outputfile):
        os.remove(outputfile)
    # obtain kmer size from hash table
    kmer_size = get_kmer_size(hash_table)
    # loop through proteins and predict structure
    for protein_index, input_seq_name in enumerate(protein_names):
        input_seq = proteins[protein_index]

        # fix input - add glycine where * occurs
        input_seq = clean_input_sequences(input_seq)
        
        # instantiate the object and run the algorithm
        if args.use_clusters:
            new_obj = DebruijnExtend(kmer_clust)
        else:
            new_obj = DebruijnExtend()
        print(f"input: \n{input_seq},\n k_mer size: {kmer_size}")
        secondary, prob = new_obj.debruijnextend(input_seq, kmer_size, hash_table)[0]
        print(f"k={kmer_size}: {secondary}")

        # save the output
        outfile = open(outputfile, "a")
        outfile.write(f"{input_seq_name}\n{prob}\n{secondary}\n")

if __name__ == "__main__":
    main()
