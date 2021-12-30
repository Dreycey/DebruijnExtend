"""
This module contains methods/objects that facilitate
basic operations.
"""
# std pkgs
import numpy as np
import random
from typing import Dict, List, Optional, Union
from pathlib import Path
import pickle
# non-std pkgs
import matplotlib.pyplot as plt




def clean_input_sequences(input_seq):
    """
    This method cleans all input sequences to ensure they will
    be compatible with the precomputed hash table. 
    """
    seq_list = []
    for aa in input_seq:
        if aa not in ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]:
            print(aa)
       
        if aa == "*":
            amino_chosen = "G"
        elif aa == "B":
            amino_chosen = np.random.choice(["N", "D"], 1, p=[0.5, 0.5])[0]
        elif aa == "Z":
            amino_chosen = np.random.choice(["Q", "E"], 1, p=[0.5, 0.5])[0]
        elif aa == "J":
            amino_chosen = np.random.choice(["L", "I"], 1, p=[0.5, 0.5])[0]
        elif aa == "X":
            amino_chosen = random.choice(["A", "R", "N", "D", "C", 
                                          "Q", "E", "G", "H", "I", 
                                          "L", "K", "M", "F", "P", 
                                          "S", "T", "W", "Y", "V"])[0]
        else:
            amino_chosen = aa
        seq_list.append(amino_chosen)
    return ''.join(seq_list) #+ input_seq[kmer_size+1:]

def readFasta(fasta_file_path: Union[str, Path]):
    """
    This function reads a fasta file

    Parameters
    ----------
    fasta file path: string OR Path

    Returns
    -------
    proteins : array of protein sequence (ordered)
    protein_names : array of protein names (ordered)
    """
    proteins, protein_names = [], []
    with open(fasta_file_path) as fasta_file:
        fasta_file_array = fasta_file.readlines()
        for line_count, fasta_line in enumerate(fasta_file_array):
            if (fasta_line[0] == ">"):
                name = fasta_line.strip("\n")
                protein_names.append(name)
                proteins.append(protein_seq) if line_count > 0 else None
                protein_seq = "" # renew sequence everytime fasta name is added.
            else: 
                protein_seq += fasta_line.strip("\n")
        proteins.append(protein_seq)
    return proteins, protein_names

def get_kmer_size(hash_table) -> int:
    """
    This function extracts the kmer size from 
    the hash table.
    """
    kmer_size = 0
    with open(hash_table, "rb") as hash_tb:
        hash = pickle.load(hash_tb)
    kmer_size = len(list(hash.keys())[0])
    return kmer_size
