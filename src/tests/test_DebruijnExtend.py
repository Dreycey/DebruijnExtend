#! usr/bin/python3
"""
This test module tests the DebruijnExtend module of DebruijnExtend
"""
# std-packages
import unittest
from pathlib import Path
from typing import Dict
# non-std library 
import numpy as np
# import relative in-house APIs
#from .. import debruijnextend
from src.debruijnextend.KmerCluster import KmerCluster
from src.debruijnextend.DebruijnExtend import DebruijnExtend, HashTableType
from src.debruijnextend.utils import readFasta, get_kmer_size, clean_input_sequences




# globals
cluster_pickle: Path = './cluster_file.pickle'
hash_table_path: Path = './prot_hashtables/prothashtable_10.p'
test_seq_1: Path = './examples/gfp.fasta'
truth = 'CCCHHHHHCCEEEEEEEEEEEECCEEEEEEEEEEEECCCCEEEEEEEECCCCCCCCHHHHCCCCCCHHHCECCHHHHHHCHHHHCCCCCEEEEEEEEECCCCEEEEEEEEEEECCEEEEEEEEEEECCCCCCCCCCCCECCCCCCEEEEEEEEHHHCEEEEEEEEEEEECCCCEEEEEEEEEEEECCCCCCCCCCCEEEEEEEEEECCCCCCCCEEEEEEEEEEECCCCCCCCCCCCC'
HashTableType = Dict[str, Dict[str, float]]
# tests
class Test_DebruijnExtend(unittest.TestCase):
    """ These tests test the DebruijnExtend algorithm"""

    def setUp(self):
        print("SetUp() method")
        kmer_clust = KmerCluster.init_struct(hash_table_path, cluster_pickle)
        self.debext = DebruijnExtend()
        self.debext_clusters = DebruijnExtend(kmer_clust)
        # obtain kmer size from hash table
        self.kmer_size = get_kmer_size(hash_table_path)

    def tearDown(self) -> None:
    #     pass
        return 0
        
    def test_get_fist_struct(self):
        """
        Test Def: 
        ---------
            Tests the get_fist_struct() method to ensure
            it reuturns a structure even if unkown.

        Definitions:
        ------------
            Sequence = 'ASKGEELFTGVVPIL'
            kmers (k=10) = ['ASKGEELFTG', 'SKGEELFTGV', 'KGEELFTGVV', 'GEELFTGVVP', 'EELFTGVVPI']

        Expected Functionality:
        ----------------------
            If a kmer is not found, then search for next possible kmer and prepend
            all ss3 possibilities. (i.e. starting points give all permutations)
        """
        k_size = 10
        sequence = 'ASKGEELFTGVVPIL'
        # get data structures ready
        prot2secondary_kmers = {'SKGEELFTGV' : {'HHHHHHCCCC' : 1}}
        primary_seq_kmers = [sequence[i:i+k_size] for i in range(len(sequence)-k_size)]
        for kmer in primary_seq_kmers:
            if kmer not in prot2secondary_kmers: prot2secondary_kmers[kmer] = {}
        # test
        possible_starts = self.debext.get_first_struct(prot2secondary_kmers, primary_seq_kmers)
        expected = [i+'HHHHHHCCC' for i in ['E', 'C', 'H']] # must be length 10
        self.assertCountEqual(list(possible_starts.keys()), expected)

    def test_accuracy(self):
        """
        Test Def: 
        ---------
            Integration test to ensure the accuracy is above the expected
            threshold for the test sequences.  

        Definitions:
        ------------

        Expected Functionality:
        ----------------------
            For a set of known proteins, the goal is to report the accuracies
            as being above the threshold for accuracy.
        """
        proteins, protein_names = readFasta(test_seq_1)
        input_seq = clean_input_sequences(proteins[0])
        secondary, prob = self.debext.debruijnextend(input_seq, 
                                                     self.kmer_size, 
                                                     hash_table_path)[0]
        print(f"k={self.kmer_size}: {secondary}")

        accuracy, conf_matrix = self.prediction_rank(secondary, truth)

        self.assertGreater(accuracy, 0.85, "The accuracy was not up to par.")

    def prediction_rank(self, prediction, true):                                 
        """                                                                         
        This function scores the accuracy of a prediction.     
                                                                                    
        INPUT: 1. an input prediction sequence
            2. the true secondary structure sequence
        OUTPUT: 1. percent guessed correctly (range: 0-1)
                2. confusion matrix
        """
        # init
        confusion_matrix = [[0,0,0],[0,0,0],[0,0,0]]
        map_dict = {"C" : 0, "E" : 1, "H" : 2}
        prediction_array = np.zeros((len(true)))
        # iterate through predicted and true
        for index in range(len(prediction)):
            # update prediction_array
            if prediction[index] == true[index]:
                prediction_array[index] = 1
            # update confusion matrix
            if (true[index] in map_dict.keys()) and (prediction[index] in map_dict.keys()):
                true_ind = map_dict[true[index]]
                pred_ind = map_dict[prediction[index]]
                confusion_matrix[true_ind][pred_ind] += 1
        # calculate the accuracy
        accuracy = round(sum(prediction_array) / len(prediction_array), 4)

        return accuracy, confusion_matrix 