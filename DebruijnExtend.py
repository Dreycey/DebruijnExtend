#! usr/bin/python3
"""
This is the primary script and code for running the DebruijnExtend algorithm
to predict the three-state secondary structure of an input fasta file with protein entries. 
DebruijnExtend can be ran in both macOS and most Linux distributions. 
For more information, please visit: https://github.com/Dreycey/DebruijnExtend
"""
# std pkgs
import math
import sys
from typing import Dict, List, Optional, Union
import numpy as np
from pathlib import Path
import os
import pickle
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
from src.debruijnextend.csvtohash import ProteinHash



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
    parser.add_argument("-ht", "--hash_table", help="the input hash table (python pickle file) \
                                                     OR kmer (int)", required=True)
    parser.add_argument("-o", "--output_file", help="the output file prefix", required=True)
    parser.add_argument("-t", "--threads", help="number of threads", required=False)
    parser.add_argument("-c", "--use_clusters", help="Path to cluster pickle OR cutoff threshold (int), \
                                                      O(M*N) full search if not used", required=False)
    parser.add_argument("-d", "--csv_data_path", help="path to CSV file containing amino-ss3 mapping (used with -ht)", required=False)
    return parser.parse_args(argv)

# RUN SCRIPT
def main():
    """
    This function controls the flow of the script.
    """
    args = parseArgs(sys.argv[1:])
    # get input ready
    proteins, protein_names = readFasta(args.input)
    outputfile = Path(args.output_file)
    current_dir_path = Path(os.path.realpath(__file__)).parent
    # updating CSV data file to use:
    if args.csv_data_path:
        data_path = Path(args.csv_data_path)
    else:
        data_path = current_dir_path / Path('data/TRAINING.csv')
    # create default hash if it doesn't exists, report error if path doesn't exist
    if args.hash_table.isdigit():
        # setting up paths and default hash table file naming
        hash_directory = current_dir_path / Path('HashtableData/')
        k_size = int(args.hash_table)
        new_hashtable_pickle = hash_directory / Path(f"hashtable_k{k_size}.pickle")
        # create hash table
        if not os.path.exists(hash_directory): os.makedirs(hash_directory)
        if not os.path.exists(new_hashtable_pickle):
            with open(new_hashtable_pickle,'wb') as outfile_hashtable:
                prothashOBJ = ProteinHash(data_path, k_size)
                prothashtable = prothashOBJ.construct_hash()
                pickle.dump(prothashtable, outfile_hashtable)
                del prothashOBJ
        hash_table = Path(new_hashtable_pickle)
    else:
        hash_table = Path(args.hash_table)
    # obtain kmer size from hash table
    kmer_size = get_kmer_size(hash_table)
    # get protein clusters
    if args.use_clusters:
        if args.use_clusters.isdigit():
            # setting up paths and default cluster file naming
            cluster_directory = current_dir_path / Path('ClusterData/')
            cut_off = int(args.use_clusters)
            new_cluster_pickle = cluster_directory / Path(f"cluster_k{kmer_size}_CT{cut_off}.pickle")
            # creating or using default names cluster file
            if not os.path.exists(cluster_directory): os.makedirs(cluster_directory)
            kmer_clust = KmerCluster.init_struct(hash_table, new_cluster_pickle, cut_off)
        else:
            kmer_clust = KmerCluster.init_struct(hash_table, args.use_clusters)
    # delete output file if exists
    if os.path.exists(outputfile):
        os.remove(outputfile)
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
