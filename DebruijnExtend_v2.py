#! usr/bin/python3
# std pkgs
import math
import pickle
import sys
# non-std pkgs
import matplotlib.pyplot as plt
from tqdm import tqdm

"""
This is the primary script and code for running the DebruijnExtend algorithm.

USAGE
         python DebruijnExtend_v2.py <input fasta> <kmer size> <output file>
EXAMPLE:
         python DebruijnExtend_v2.py examples/gfp.fasta 4 gfp.ss3

Updates (v2; 11/28/2020): 
    * Here unelongated sequences are not extended further.
    * There is a max size set for the debruijn extension. This helps
      speed the algorithm dramatically and allows for scaling to big
      protein sequences.
"""



class DebruijnExtend():
    """
    This class impliments a debruijn extend 
    """
    
    # CONSTRUCTOR
    def __init__(self):
        """
        Constructor

        This will initialize the sequence of the protein, as well as the 
        predicted 3-state secondary structure. 
        """
        self.seq_name = ""
        self.sequence = ""
        self.secoondary = ""

    # get kmers 
    def get_kmers(self, primary_prot_seq, kmer_length):
        """
        returns an array of the input string as chunks length k
        """
        kmer_array = [primary_prot_seq[i:i+kmer_length] 
                      for i in range(len(primary_prot_seq)-(kmer_length-1))]  

        return kmer_array

    def find_secondary_structs(self, primary_karray, hash_table):
        """
        This method returns only the relevant hashes from the larger hash table.
        """
        potential_secondaries = {} # this will have all layers for our graph
        for kmer in primary_karray:
            # turn into probalities 
            temp_dict = {}
            summ = 0
            for key, value in hash_table[kmer].items():
                summ += value
                temp_dict[key] = value
            for key, value in temp_dict.items():
                temp_dict[key] = (-1) * math.log(round(value / summ, 5))    
            # add potential secondarys structure to new dictionary
            potential_secondaries[kmer] = temp_dict

        return potential_secondaries

    def stitchextend(self, primary_karray, potential_secondaries, k):
        """
        CORE ALGORITHM.

        This method is the core of the debruijnextend algorithm. It impliments a
        dynamic programming approach to traverse through each layer of the 
        debruijn graph, then uses edge contraction to create the new input to 
        the next iteration.
        """
        # Initialize with first
        stitchextend_dict = potential_secondaries[primary_karray[0]].copy()
        # make max size constraint for stitchextend_dict
        max_dict_size = 10000
 
        print("starting 1") 
        # LOOP 1: looping through the layers
        for kmer_i in tqdm(range(1,len(primary_karray))):
            kmer_lay = primary_karray[kmer_i]
            seq_ext = {}
            stitchextend_dict_iplus = {}
            # LOOP 2: loop through kmers per layer
            for kmer_in_layer, kmer_prob in potential_secondaries[kmer_lay].items():
                # LOOP 3: loop through extended sequences
                print(f"number of sequences: {len(stitchextend_dict.keys())}")
                for seq, prob in stitchextend_dict.items():
                    end_of_ext = seq[-(k-1):] # last k-1
                    start_of_kmer = kmer_in_layer[:k-1] # up to k-1
                    seq_ext[seq] = 0 # delete on next layer #TODO: store in
                                     # final if not extended.
                    if end_of_ext == start_of_kmer:

                        # if equal, THEN STITCH AND EXTEND
                        extended_seq = seq + kmer_in_layer[-1]
                        extended_prob = prob + kmer_prob
                        # add new seq and prob
                        stitchextend_dict_iplus[extended_seq] = extended_prob

            print(f"number of sequences before: {len(stitchextend_dict.keys())}")
            for sequence_extended in seq_ext.keys():
                del stitchextend_dict[sequence_extended] # delete the old
            print(f"number of sequences after del: {len(stitchextend_dict.keys())}")
            top_probable_seqs_list = sorted(stitchextend_dict_iplus.items(), 
                                       key=lambda item: item[1])[:max_dict_size]
            top_probable_seqs_dict = {k: v for k, v in top_probable_seqs_list}
            stitchextend_dict.update(top_probable_seqs_dict) # add the new
            print(f"number of sequences after add: {len(stitchextend_dict.keys())}")

        return stitchextend_dict

    # Algorithm
    def debruijnextend_v1(self, primary_seq, k):
        """
        This method takes in a primary protein sequence annd returns the
        secondary structure predicted with the highest probability.
        INPUT: primary protein sequence (type: string), 
               pickled hash table (type: python dict),
               kmer length (type: int)
        OUTPUT: 3-based secondary structure (using CEH, type: string) 
        """
        primary_karray = self.get_kmers(primary_seq, k)
        hash_table = pickle.load(open(f"prot_hashtables/prothashtable_{k}.p", "rb"))

        ###
        # STEP 1: Find corresponding secondary structures
        ###
        print(f"STEP 1: Hashing Method. Find corresponding secondary structures: \n")
        potential_secondaries = self.find_secondary_structs(primary_karray,
                                                            hash_table)

        ###
        # STEP 2: Connect the layers, use dynamic programming per layer (i.e.
        # StitchExtend)
        ###
        print(f"STEP 2: StitchExtend Method. Looping through layers: \n")
        stitchextend_dict = self.stitchextend(primary_karray,
                                              potential_secondaries, k)

        ####
        # STEP 3: choose the top 10
        ####
        print(f"STEP 3: Choosing top 50 predictions \n")
        top_number = 50
        out_array = sorted(stitchextend_dict.items(), 
                           key=lambda item: item[1])[:top_number]

        return out_array

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
    k = int(sys.argv[2])
    outputfile = sys.argv[3]
    
    # instantiate the object
    new_obj = DebruijnExtend()

    # run the algorithm
    print(f"input: \n {input_seq},\n k_mer size: {k}")
    secondary, prob = new_obj.debruijnextend_v1(input_seq, k)[0]
    print(f"k={k}: {secondary}")

    # save the output
    outfile = open(outputfile, "w")
    outfile.write(f"{input_seq_name} \n{prob} \n{secondary}")

if __name__ == "__main__":
    main()
