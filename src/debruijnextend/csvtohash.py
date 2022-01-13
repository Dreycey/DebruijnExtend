#! usr/bin/python3
"""
This module creates a hash table with sequences and the corresponding
secondary structures.

USAGE:
    python csvtohash.py secondarystructure_COMP533.csv prothashtable.p

OUTPUT:
    GENERAL OUTPUT:
        hashtable = { "sequence" : {"secondary structure" : counts} }   

    SPECIFIC OUTPUT EXAMPLE:
        hashtable = {"AGCTH" : {"CEEEH" : 2, ..., "EEEEE" : 30},
                      ...
                     "TTPGH" : {"HHHEE" : 34, ..., "HHEEE" : 2}}
"""
from typing import Dict
import sys
import pandas as pd
from tqdm import tqdm
import pickle
from pathlib import Path



HashTableType = Dict[str, Dict[str, float]]
# class for the hash table
class ProteinHash:
    """
    This object holds the mappings from sequence to potential 
    secondary structure. This allows for a fast mapping of 
    sequence to structure for building a graph.
    """
    
    def __init__(self, csv_in, kmer_size):
        self.csv_file: Path = Path(csv_in)
        self.prothash = None # must construct on initialization
        self.k_mersize = kmer_size 
        self.sequence_column = 2
        self.structure_column = 3
        
    def parse_csv(self, csv_file=None, delimitor: str=",") -> pd.DataFrame:
        """
        This method parses an input csv and places the sequences
        and corresponding secondary structures into a pandas dataframe.

        Parameters
        ----------
        delimitor : str [DEFAULT: ","]
            The delimitor to use per row of the input file.

        Returns
        -------
        potential_secondaries: pd.DataFrame
            Returns a dataframe connecting primary sequences with the 
            corresponding secondary structure.
        """
        if csv_file == None:
            csv_file = self.csv_file
        with open(csv_file) as opened_file:
            file_by_lines = opened_file.readlines()
            # construct the sequence and structure vectors from csv
            sequence_vector = []
            structure_vector = []
            for line_number, csv_row in enumerate(file_by_lines):
                if line_number > 0:
                    csv_row_list = csv_row.strip("\n").split(delimitor)
                    sequence_vector.append(csv_row_list[self.sequence_column].strip("'"))
                    structure_vector.append(csv_row_list[self.structure_column].strip("'"))
            # add sequences and structures to a pandas dataframe
            df = {'sequence': sequence_vector, 'structure': structure_vector}
            structseq_dataframe = pd.DataFrame(data=df)

        return structseq_dataframe
    
    def construct_hash(self, seq_column_name: str="sequence", 
                             struct_column_name: str="structure") -> HashTableType:
        """
        This method creates the hash table used by the DebruijnExtend algorithm.

        Parameters
        ----------
        seq_column_name : str [DEFAULT: "sequence"]
            The column name for the sequences in a pandas dataframe
        struct_column_name : str [DEFAULT: "structure"]
            The column name for the structures in a pandas dataframe

        Returns
        -------
        potential_secondaries: HashTableType
            returns a hashable datastructure.
        """
        prothash = {}
        protein_df = self.parse_csv()
        print(f"Creating hash table for kmer size = {self.k_mersize} (takes a couple minutes)")
        # make the hash table
        for kmer_start_position, row in tqdm(protein_df.iterrows()):
            primary_sequence, secondary_structure = row[seq_column_name], row[struct_column_name]
            for kmer_start_position in range(len(primary_sequence) - (self.k_mersize-1)):
                kmer_end_position = kmer_start_position + self.k_mersize
                protein_seq_at_index = primary_sequence[kmer_start_position:kmer_end_position]
                structure_seq_at_index = secondary_structure[kmer_start_position:kmer_end_position]
                prothash = self.update_protein_hash_table(prothash,
                                                          protein_seq_at_index,
                                                          structure_seq_at_index)
        return prothash
    
    def update_protein_hash_table(self, prothash: HashTableType, 
                                        protein_seq_at_index:str, 
                                        structure_seq_at_index:str) -> HashTableType:
        """
        This method updates the protein hash table.
        """
        if "*" not in protein_seq_at_index: # TODO: check what asterisk is!!!!!!!
            if protein_seq_at_index in prothash.keys():
                if structure_seq_at_index in prothash[protein_seq_at_index].keys():
                    prothash[protein_seq_at_index][structure_seq_at_index] += 1
                else:
                    prothash[protein_seq_at_index][structure_seq_at_index] = 1
            else:
                prothash[protein_seq_at_index] = {}
                prothash[protein_seq_at_index][structure_seq_at_index] = 1   
        return  prothash

    def map(self, input_key: str) -> Dict[str, float]:
        """
        return the values for the input key.

        Parameters
        ----------
        input_key: str
            a sequence of length k, where k is the sliding window length

        Returns
        -------
        the potential secondary structures and associated counts.   
        """
        return self.prothash[input_key]

def main():
    # defining the input arguments
    csvfile = "/Users/dreyceyalbin/Dropbox/Fall2020-classes/Algorithms/project/DebruijnExtend/py_scripts/training_1.csv"; #sys.argv[1]
    outpickle_name = "out" #sys.argv[2] 
    kmer_size = 4 # int(sys.argv[3])
    # construct the hash table
    prothashOBJ = ProteinHash(csvfile, kmer_size)
    prothashtable = prothashOBJ.construct_hash()
    print(prothashtable)
    # pickle the hash table
    outfile = open(outpickle_name,'wb')
    pickle.dump(prothashtable, outfile)

if __name__ == "__main__":
    main()
