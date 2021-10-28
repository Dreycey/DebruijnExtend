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
import pickle
import sys
import os
from typing import Dict, List, Union
import numpy as np
# non-std pkgs
import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path
import random


HashTableType = Dict[str, Dict[str, float]]

class DebruijnExtend():
    """
    This class impliments the debruijn extend algorithm.
    """
    
    def __init__(self):
        """
        Constructor

        This will initialize the sequence of the protein, as well as the 
        predicted 3-state secondary structure. 
        """
        self.seq_name = ""
        self.sequence = ""
        self.secondary = ""
        self.secondary_alphabet = ["C","E","H"]

    def get_kmers(self, primary_prot_seq: str, kmer_length: int) -> List[str]:
        """
        Returns an array of the input string as chunks length k.

        Parameters
        ----------
        primary_prot_seq : str
            The primary seqeunce for the protein.
        kmer_length : int
            The kmer length for the sliding window.

        Returns
        -------
        kmer_array: list
            The output kmers for a given input string.
        """
        kmer_array = [primary_prot_seq[i:i+kmer_length] 
                      for i in range(len(primary_prot_seq)-(kmer_length-1))]  

        return kmer_array

    def find_secondary_structs(self, primary_karray: list, 
                                     hash_table: HashTableType) -> HashTableType:
        """
        This method returns only the relevant hashes from the larger hash table. 

        Parameters
        ----------
        primary_karray : list
            The primary seqeunce for the protein.
        hash_table: HashTableType
            The INPUT hash table with mappings from kmer to secondary 
            structure kmer and observed COUNTS.

        Returns
        -------
        potential_secondaries: HashTableType
            The OUTPUT hash table with mappings from kmer to secondary 
            structure kmer and observed negative log PROBABILITIES.
        """
        potential_secondaries = {} # this will have all layers for our graph
        for kmer in primary_karray:
            temp_dict = {}
            if kmer in hash_table.keys(): 
                for secondary_kmer, observed_count in hash_table[kmer].items():
                    summ = np.sum([observed_count for observed_count in hash_table[kmer].values()])
                    temp_dict[secondary_kmer] = (-1) * math.log(round(observed_count / summ, 20)) # turn into probalities 
            potential_secondaries[kmer] = temp_dict 

        return potential_secondaries

    def stitchextend(self, primary_seq_kmers: List[str], 
                           primary2secondary_kmers: HashTableType, 
                           kmer_size: int = 4, 
                           max_dict_size: int = 100,
                           prob_cutoff: float = .50) -> Dict[str, float]:
        """
        The stitch-extend method is the core algorithm in DebruijnExtend. 

        Description
        ----------
        This method is the core of the debruijnextend algorithm. It impliments a
        dynamic programming approach to traverse through each layer of the 
        debruijn graph, then uses edge contraction to create the new input to 
        the next iteration.

        Parameters
        ----------
        primary_seq_kmers : list
            The list of kmers in the primary sequences (needs to be the same size as the kmers in the 
            primary2secondary_kmers table)
        primary2secondary_kmers: HashTableType
            The hash table with mappings from primary kmer to secondary 
            structure kmer and observed negative log PROBABILITIES.
        kmer_size: int [DEFAULT: 4]
            This is the kmer size used for the algorithm. This is will have corresponding 
            table kmer size (in primary2secondary_kmers) and a kmer list of the primary 
            sequence with the correct size (in primary_seq_kmers).
        max_dict_size: int [DEFAULT: 100]
            This is the max dictionary size, which indicates the max number of extended
            sequences evaluated per iteration. Setting a size limit for this paramater
            prevents the time complexity from blowing up.
        prob_cutoff: float [DEFAULT: 0.10]
            The probabilitiy cutoff ensures that there are sequences for each iteration. If no
            extended sequences are found, then the heuristic is to add all possible combintations
            to the possible outputs.

        Returns
        -------
        stitchextend_dict: Dict[str, float]
            The final output is a dictionary of potential secondary structure sequences 
            with cognate negative log proabilities. 
        """
        # Initialize with first
        stitchextend_dict = primary2secondary_kmers[primary_seq_kmers[0]].copy()
        # LOOP 1: looping through the layers
        for kmer_i in tqdm(range(1,len(primary_seq_kmers))):
            protein_kmer_at_current_layer = primary_seq_kmers[kmer_i]
            stitchextend_dict_iplus = {}
            # LOOP 2: loop through kmers per layer
            for secondary_kmer_in_layer, secondary_kmer_probability in \
                                  primary2secondary_kmers[protein_kmer_at_current_layer].items():
                # LOOP 3: loop through extended sequences
                for extended_sequence, extended_probability in stitchextend_dict.items():
                    end_of_extended_sequence = extended_sequence[-(kmer_size-1):] # last k-1
                    start_of_kmer = secondary_kmer_in_layer[:-1] # up to k-1
                    if end_of_extended_sequence == start_of_kmer: # if equal, THEN STITCH AND EXTEND
                        extended_seq = extended_sequence + secondary_kmer_in_layer[-1]
                        extended_prob = extended_probability + secondary_kmer_probability
                        stitchextend_dict_iplus[extended_seq] = extended_prob
            # LOOP 2 END: Update the stitch-extend dictionary
            stitchextend_dict = self.apply_heurstics(stitchextend_dict, 
                                                     stitchextend_dict_iplus, 
                                                     prob_cutoff, 
                                                     max_dict_size)
        return stitchextend_dict

    def apply_heurstics(self, stitchextend_dict,
                              stitchextend_dict_iplus, 
                              prob_cutoff, 
                              max_dict_size):
        """
        This method applies emperically-derived heurstics to the DebruijnExtend algorithm.

        Description
        ----------
        Heuristic 1:
            If prob_cutoff*100 % of max_dict_size not added, then extend
            previous with all 3 C,E,H with probability 1.   
        Heuristic 2: 
            Only keep the extended sequences with lowest neg log probability (highest probability).
        """
        # Heuristic 1
        if len(stitchextend_dict_iplus.keys()) < prob_cutoff*max_dict_size:
            stitchextend_dict_iplus = {} # empty the hash
            for extended_sequence, extended_probability in stitchextend_dict.items():
                for nucleo_ext in self.secondary_alphabet:
                    extended_seq = extended_sequence + nucleo_ext
                    stitchextend_dict_iplus[extended_seq] = extended_probability # adding zero, since -log(1)=0 
        # Heuristic 2
        top_probable_seqs_list = sorted(stitchextend_dict_iplus.items(), 
                                    key=lambda item: item[1])[:max_dict_size]
        return {k: v for k, v in top_probable_seqs_list}

    def debruijnextend(self, primary_seq: str, kmer_size: int) -> List[str]:
        """
        This method takes in a primary protein sequence annd returns the
        secondary structure predicted with the highest probability.

        Parameters
        ----------
        primary protein sequence: string
        kmer length: int

        Returns
        -------
        3-based secondary structure (using CEH, type: string) 
        """
        primary_karray = self.get_kmers(primary_seq, kmer_size)
        curr_loc = os.path.abspath(os.path.dirname(__file__))
        hash_table = pickle.load(open(f"{curr_loc}/prot_hashtables/prothashtable_{kmer_size}.p", "rb"))

        print(f"STEP 1: Hashing Method. Find corresponding secondary structures: \n")
        potential_secondaries = self.find_secondary_structs(primary_karray, hash_table)

        print(f"STEP 2: StitchExtend Method. Looping through layers /Dynamic Programming: \n")
        stitchextend_dict = self.stitchextend(primary_karray, potential_secondaries, kmer_size)

        print(f"STEP 3: Choosing top 50 predictions \n")
        top_number = 50
        out_array = sorted(stitchextend_dict.items(), key=lambda item: item[1])[:top_number]

        return out_array

def clean_input_sequences(input_seq, kmer_size):
    """
    This method cleans all input sequences to ensure they will
    be compatible with the precomputed hash table. 
    """
    seq_list = []
    for aa in input_seq[:kmer_size+1]:
        if aa == "*":
            seq_list.append("G")
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
    return ''.join(seq_list) + input_seq[kmer_size+1:]

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


###
# RUN SCRIPT
###
def main():
    """
    This function controls the flow of the script.
    """
    # get input ready
    proteins, protein_names = readFasta(sys.argv[1])
    kmer_size = int(sys.argv[2])
    outputfile = sys.argv[3]
    # delete output file if exists
    if os.path.exists(outputfile):
        os.remove(outputfile)
    # loop through proteins and predict structure
    for protein_index, input_seq_name in enumerate(protein_names):
        input_seq = proteins[protein_index]

        # fix input - add glycine where * occurs
        input_seq = clean_input_sequences(input_seq, kmer_size)
        
        # instantiate the object and run the algorithm
        new_obj = DebruijnExtend()
        print(f"input: \n {input_seq},\n k_mer size: {kmer_size}")
        secondary, prob = new_obj.debruijnextend(input_seq, kmer_size)[0]
        print(f"k={kmer_size}: {secondary}")

        # save the output
        outfile = open(outputfile, "a")
        outfile.write(f"{input_seq_name} \n {prob} \n {secondary} \n")

if __name__ == "__main__":
    main()
