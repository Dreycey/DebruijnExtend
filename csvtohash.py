#! usr/bin/python3
import numpy
import sys
import pandas as pd
from tqdm import tqdm
import pickle

"""
This module creates a hash table with sequences and the corresponding
secondary structures. Example:

USAGE:
python csvtohash.py secondarystructure_COMP533.csv prothashtable.p

OUTPUT:

                           GENERAL OUTPUT:
hashtable = {                                                                   
             "sequence" : {                                                        
                           "secondary structure" : counts,                           
                          }
            }   

                           SPECIFIC OUTPUT EXAMPLE:
hashtable = {
             "AGCTH" : {
                        "CEEEH" : 2, 
                        "EEEEE" : 30,  
                       }
             "TTPGH" : {                                                        
                        "HHHEE" : 34,                                            
                        "HHEEE" : 2,                                           
                       } 
             "WGGTT" : {                                                        
                        "HHHHH" : 50,                                            
                        "EEHHH" : 4,                                           
                       }  
            }

"""

# class for the hash table
class proteinhash:
    """
    This object holds the mappings from sequence to potential 
    secondary structure. This allows for a fast mapping of 
    sequence to structure for building a graph.
    """
    
    def __init__(self, csv_in, k):
        # define constructor variables.
        self.csv_file = open(csv_in)
        self.prothash = None # must construct on initialization
        self.k_mersize = k # sliding window size for data structure

        # construct the hash table
        self.construct_hash()
        
    def parse_csv(self):
        """
        This method parses an input csv and send the information into 
        a seperate method which constructs the hash table.
        """
        file_by_lines = self.csv_file.readlines()
        # construct the sequence and structure vectors from csv
        sequence_vector = []
        structure_vector = []
        for line in file_by_lines:
            line = line.split(",")
            sequence_vector.append(line[2])
            structure_vector.append(line[4])
        # add sequences and structures to a pandas dataframe
        df = {'sequence': sequence_vector, 'structure': structure_vector}
        structseq_dataframe = pd.DataFrame(data=df)
        
        return structseq_dataframe
    
    def construct_hash(self):
        """
        This meth creates a hash table.

        OUTPUT: returns a hashable datastructure.
        """
        prothash = {}
        protein_df = self.parse_csv() # creates df needed for building hash table
        k = self.k_mersize
        
        # make the hash table
        for index, row in tqdm(protein_df.iterrows()):
            seq = row['sequence']
            struct = row['structure']
            for index in range(len(seq) - (k-1)):
                start_pos = index
                end_pos = index + k
                seq_sub = seq[start_pos:end_pos]
                struct_sub = struct[start_pos:end_pos]
                if "*" not in seq_sub: # TODO: check what asterisk is!!!!!!!
                    if seq_sub in prothash.keys():
                        if struct_sub in prothash[seq_sub].keys():
                            prothash[seq_sub][struct_sub] += 1
                        else:
                            prothash[seq_sub][struct_sub] = 1
                    else:
                        prothash[seq_sub] = {}
                        prothash[seq_sub][struct_sub] = 1 
            
        return prothash

    def map(self, input_key):
        """
        return the values for the input key.

        INPUT: a sequence of length k, where k is the sliding window length
        """
        return self.prothash[input_key]



# main for testing the hash object.
def main():
    # defining the input arguments
    csvfile = sys.argv[1]
    outpickle_name = sys.argv[2] 
    kmer_size = int(sys.argv[3])
    # construct the hash table
    prothashOBJ = proteinhash(csvfile, kmer_size)
    prothashtable = prothashOBJ.construct_hash()
    # pickle the hash table
    outfile = open(outpickle_name,'wb')
    pickle.dump(prothashtable,outfile)

if __name__ == "__main__":
    main();
