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

Updates (10/04/2021): 
    * refactoring all methods/functions
        * Adding type hinting, imrpoved doc strings
        * getting rid of other "versions"
"""
# std pkgs
import math
import pickle
import sys
import os
from typing import Dict, List
import numpy as np
import multiprocessing as mp
# non-std pkgs
import matplotlib.pyplot as plt
from tqdm import tqdm




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

    def find_secondary_structs(self, primary_karray: list, hash_table: HashTableType) -> HashTableType:
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
            else:
                temp_dict = self.get_close_kmers(hash_table, kmer)
            potential_secondaries[kmer] = temp_dict 

        # pool = mp.Pool(6)
        # pool.starmap(self.find_secondary, [(kmer, potential_secondaries, hash_table) for kmer in tqdm(primary_karray)])
        # pool.join()
        return potential_secondaries

    def find_secondary(self, prot_kmer, potential_secondaries, hash_table):
        """
        given kmer, find secondary structures. return dict.
        """
        temp_dict = {}
        if prot_kmer in hash_table.keys(): 
            for secondary_kmer, observed_count in hash_table[prot_kmer].items():
                summ = np.sum([observed_count for observed_count in hash_table[prot_kmer].values()])
                temp_dict[secondary_kmer] = (-1) * math.log(round(observed_count / summ, 20)) # turn into probalities 
        else:
            temp_dict = self.get_close_kmers(hash_table, prot_kmer)
        potential_secondaries[prot_kmer] = temp_dict 

    def stitchextend(self, primary_seq_kmers: List[str], 
                           prot2secondary: HashTableType, 
                           kmer_size: int = 4, 
                           max_dict_size: int = 100,
                           prob_cutoff: float = .10) -> Dict[str, float]:
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
        stitchextend_dict = prot2secondary[primary_seq_kmers[0]].copy()

        # LOOP 1: looping through the layers
        for kmer_i in tqdm(range(1,len(primary_seq_kmers))):
            protein_kmer = primary_seq_kmers[kmer_i]
            stitchextend_dict_iplus = {}
            print(f"stitchextend_dict =  {stitchextend_dict}")
            # LOOP 2: loop through kmers per layer 
            if protein_kmer in prot2secondary:
                print("protein kmer found!")
                secondary_kmer = prot2secondary[protein_kmer]
                print(f"OUTPUT: {secondary_kmer}")
            else:
                print("finding close kmers")
                secondary_kmer = self.get_close_kmers(prot2secondary, protein_kmer)
            print(f"secoondaary structures: {secondary_kmer}")
            for secondary_k, secondary_prob in secondary_kmer.items():
                # LOOP 3: loop through extended sequences
                for extended_sequence, extended_probability in stitchextend_dict.items():
                    end_of_extended_sequence = extended_sequence[-(kmer_size-1):] # last k-1
                    start_of_kmer = secondary_k[:-1] # up to k-1
                    if end_of_extended_sequence == start_of_kmer: # if equal, THEN STITCH AND EXTEND
                        extended_seq = extended_sequence + secondary_k[-1]
                        extended_prob = extended_probability + secondary_prob
                        stitchextend_dict_iplus[extended_seq] = extended_prob
            # LOOP 2 END: Update the stitch-extend dictionary
            stitchextend_dict = self.apply_heurstics(stitchextend_dict, 
                                                     stitchextend_dict_iplus, 
                                                     prob_cutoff, 
                                                     max_dict_size)
        return stitchextend_dict

    def get_close_kmers(self, prot2secondary, protein_kmer, top_N=10):
        """
        This method finds kmers that are close and uses the 
        structures for the top N.

        Output: {"HCHCG", 0.4}
        """
        def hamming_dist(k1, k2):
            val = 0
            for ind, char in enumerate(k1):
                if char != k2[ind]: val += 1
            return val

        output_dict = {}
        priority_queue = []
        highest_score = float("inf")
        print(f"\t finding close seq to {protein_kmer}")
        for protein_seq, secondary_structs in tqdm(prot2secondary.items()):
            hamming_score = hamming_dist(protein_seq, protein_kmer)
            if hamming_score < highest_score:
                priority_queue.append((hamming_score, secondary_structs))
                priority_queue.sort(key=lambda a: a[0])
            if len(priority_queue) > top_N: priority_queue.pop(-1)
            highest_score = priority_queue[-1][0]
        print(priority_queue)

        for saved_res in priority_queue:
            output_dict.update(saved_res[1])
        return output_dict

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
        print(f"before heuristic: {stitchextend_dict_iplus}")
        # Heuristic 1
        if len(stitchextend_dict_iplus.keys()) < prob_cutoff*max_dict_size:
            print("OKAY0")
            stitchextend_dict_iplus = {} # empty the hash
            for extended_sequence, extended_probability in stitchextend_dict.items():
                # extend with all H,C,E
                print("OKAY1")
                for nucleo_ext in self.secondary_alphabet:
                    print("OKAY2")
                    extended_seq = extended_sequence + nucleo_ext
                    stitchextend_dict_iplus[extended_seq] = extended_probability # adding zero, since -log(1)=0
        print(f"after heuristic: {stitchextend_dict_iplus}")
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
        else:
            seq_list.append(aa)
    return ''.join(seq_list) + input_seq[kmer_size+1:]

###
# RUN SCRIPT
###
def main():
    """
    This function controls the flow of the script.
    """
    # get input ready
    input_seq_name = open(sys.argv[1]).readlines()[0].strip("\n")
    input_seq = open(sys.argv[1]).readlines()[1].strip("\n")
    kmer_size = int(sys.argv[2])
    outputfile = sys.argv[3]

    # fix input - add glycine where * occurs
    input_seq = clean_input_sequences(input_seq, kmer_size)
    
    # instantiate the object and run the algorithm
    new_obj = DebruijnExtend()
    print(f"input: \n {input_seq},\n k_mer size: {kmer_size}")
    secondary, prob = new_obj.debruijnextend(input_seq, kmer_size)[0]
    print(f"k={kmer_size}: {secondary}")

    # save the output
    outfile = open(outputfile, "w")
    outfile.write(f"{input_seq_name} \n{prob} \n{secondary}")

if __name__ == "__main__":
    main()
