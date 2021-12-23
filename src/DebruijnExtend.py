"""
DebruijnExtend module

This module contains all of the classes/modules
necessary for the debext implimentation. 
"""
# std pkgs
from typing import Dict, List, Optional, Union
from pathlib import Path
import os
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import itertools
import math
# non-std pkgs
import pickle
from tqdm import tqdm
# in-house pkgs




HashTableType = Dict[str, Dict[str, float]]

class DebruijnExtend():
    """
    This class impliments the debruijn extend algorithm.
    """
    
    def __init__(self, kmer_clusters=None):
        """
        Constructor

        This will initialize the sequence of the protein, as well as the 
        predicted 3-state secondary structure. 
        """
        self.seq_name = ""
        self.sequence = ""
        self.secondary = ""
        self.secondary_alphabet = ["E","H","C"]
        self.kmer_clusters = kmer_clusters

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
                    #print(secondary_kmer, observed_count)
                    summ = np.sum([observed_count for observed_count in hash_table[kmer].values()])
                    temp_dict[secondary_kmer] = (-1) * math.log(round(observed_count / summ, 20)) # turn into probalities 
            else:
                if self.kmer_clusters:
                    temp_dict = self.get_close_kmers_clusters(hash_table, kmer)
                else:
                    temp_dict = self.get_close_kmers(hash_table, kmer)
            potential_secondaries[kmer] = temp_dict 


        # self.hash_table = hash_table
        # #pool = mp.Pool(8)
        # with ThreadPoolExecutor(max_workers = 10) as executor:
        #     potential_secondaries = dict(executor.map(self.find_secondary, [kmer for kmer in primary_karray]))
        # #potential_secondaries = dict(pool.starmap(self.find_secondary, [(kmer,hash_table) for kmer in primary_karray]))
        # #pool.close()
        # #pool.join()
        # for result in potential_secondaries:
        #     print(result)
        return potential_secondaries

    def get_close_kmers_clusters(self, hash_table, kmer, top_N=10):
        """
        finds possibl structures using clusers instead of all vs all
        """
        print(f"looking at kmer: {kmer}")
        def hamming_dist(k1, k2):
            val = 0
            for ind, char in enumerate(k1):
                if char != k2[ind]: val += 1
            return val
        # find related clusters    
        print("looking at clusters")
        threshold = 8
        kmers_to_look_at = []
        for centroid, cluster_kmers in tqdm(self.kmer_clusters.clusters.items()):
            if hamming_dist(centroid, kmer) < threshold:
                kmers_to_look_at += [kmer_i for kmer_i in cluster_kmers]
        # use found kmers for further evaluation
        priority_queue = []
        highest_score = float("inf")
        for kmer_j in tqdm(kmers_to_look_at):
            hamming_score = hamming_dist(kmer_j, kmer)
            if hamming_score < highest_score:
                secondary_structs = hash_table[kmer_j]
                priority_queue.append((hamming_score, secondary_structs))
                priority_queue.sort(key=lambda a: a[0])
            if len(priority_queue) > top_N: priority_queue.pop(-1)
            highest_score = priority_queue[-1][0]
        print(priority_queue )
        # turn into output dictionary
        output_dict = {}
        for saved_res in priority_queue:
            output_dict.update(saved_res[1])
        print(output_dict)
        return output_dict

    def find_secondary(self, prot_kmer):
        """
        given kmer, find secondary structures. return dict.
        """
        temp_dict = {}
        if prot_kmer in self.hash_table.keys(): 
            for secondary_kmer, observed_count in self.hash_table[prot_kmer].items():
                summ = np.sum([observed_count for observed_count in self.hash_table[prot_kmer].values()])
                temp_dict[secondary_kmer] = (-1) * math.log(round(observed_count / summ, 20)) # turn into probalities 
        else:
            temp_dict = self.get_close_kmers(self.hash_table, prot_kmer)
        return prot_kmer, temp_dict

    def stitchextend(self, primary_seq_kmers: List[str], 
                           prot2secondary: HashTableType, 
                           kmer_size: int = 4, 
                           max_dict_size: int = 100,
                           prob_cutoff: float = .70) -> Dict[str, float]:
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

        #stitchextend_dict = prot2secondary[primary_seq_kmers[0]].copy()
        stitchextend_dict = self.get_fist_struct(prot2secondary, primary_seq_kmers)
        print(stitchextend_dict)

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

    def get_fist_struct(self, primary2secondary_kmers, primary_seq_kmers):
        """
        This method finds the fist sequence to start extending from.
        """
        first_struct = {}
        skip_count = 0
        for kmer in primary_seq_kmers:
            struct = primary2secondary_kmers[kmer].copy()
            print(f"struct: {struct}")
            print(len(struct.keys()))
            if len(struct) == 0:
                print("EMPTY")
                skip_count += 1
            else:
                # continue
                # alphabet = ["H", "E", "C"]
                # print(f"FOUND: {struct}")
                for comb in itertools.combinations(self.secondary_alphabet, skip_count):
                    prepend_seq = ''.join(comb)
                    print(f" prepend seq: {prepend_seq}")
                    for ss3, prob in struct.items():
                        new_seq = prepend_seq+ss3
                        k = len(ss3)
                        first_struct[new_seq[:k]] = prob
                break
        return first_struct

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
            #stitchextend_dict_iplus = {} # empty the hash

            for extended_sequence, extended_probability in stitchextend_dict.items():
                # extend with all H,C,E
                print("OKAY1")
                for nucleo_ext in self.secondary_alphabet:
                    print("OKAY2")
                    extended_seq = extended_sequence + nucleo_ext
                    stitchextend_dict_iplus[extended_seq] = extended_probability + 100 # adding 100 -> low prob.
        # Heuristic 2
        top_probable_seqs_list = sorted(stitchextend_dict_iplus.items(), 
                                        key=lambda item: item[1])[:max_dict_size]
        return {k: v for k, v in top_probable_seqs_list}

    def debruijnextend(self, primary_seq: str, 
                             kmer_size: int, 
                             hash_table_path: Optional[Path]) -> List[str]:
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
        hash_table = pickle.load(open(hash_table_path, "rb"))

        print(f"STEP 1: Hashing Method. Find corresponding secondary structures: \n")
        potential_secondaries = self.find_secondary_structs(primary_karray, hash_table)

        print(f"STEP 2: StitchExtend Method. Looping through layers /Dynamic Programming: \n")
        stitchextend_dict = self.stitchextend(primary_karray, potential_secondaries, kmer_size)

        print(f"STEP 3: Choosing top 50 predictions \n")
        top_number = 50
        out_array = sorted(stitchextend_dict.items(), key=lambda item: item[1])[:top_number]

        return out_array