#! usr/bin/python3
# std pkgs
import sys
import os
import numpy as np
# non-std pkgs
import pandas as pd
import pickle

"""
This script runs all benchmarking steps to evaluate the accuracy of the 
DebruijnExtend algorithm and implimentation.

USAGE:
    python benchmark_DebruijnExtend.py <benchmark config file>
EXAMPLE:
    python benchmark_DebruijnExtend.py benchmarkTEST.config
"""

def prediction_rank(prediction, true):                                 
    """                                                                         
    this function scores the accuracy of a prediction.     
                                                                                
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

def createpandas_df(csv_path):
    """
    Creates pandas dataframe using csv file.
    """
    file_by_lines = open(csv_path).readlines()
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

def parse_psipred_file(psipred_output):
    """
    This function parses the psipred output file and returns the
    secondary structure predicted.
    """
    opened_file = open(psipred_output).readlines() 
    seq = ""
    for line in opened_file:
        line = line.strip().split(" ")
        seq += line[2]
    
    return seq

def benchmark_psipredict(binary_path, test_df, test_iteration):
    """
    this function performs the benchmarking for the software: psipredict
    
    INPUT: 1. Path to the softwares binary, 
           2. test/true pandas dataframe,
           3. test iteration
    OUTPUT:
           dictionary = { 1 : [[.8, .4, ..., 0.7], # percent guessed correctly
                               [254, 223, ..., 30], # length per
                               [[], # confusion matrix
                                [],
                                []
                               ]
                              ], 
                          ...
                        }
    """
    output_dict = {test_iteration : [[], [], np.zeros((3,3))] }                                                                     
    seq_count = 0
    for index, row in test_df.iterrows(): # loop through tests
        prot_test_seq = row['sequence']
        prot_test_seq_len = len(prot_test_seq)
        prot_true_sec = row['structure']
        # input for PsiPred
        fa_filename = f"fa_filename_{test_iteration}_{seq_count}"
        fasta_name = f"fasta_name_{test_iteration}_{seq_count}"
        # run the command
        create_fasta(fa_filename, fasta_name, prot_test_seq)
        CMD = f"{binary_path} {fa_filename}"
        os.system(CMD)
        remove_file(fa_filename)
        # parse the output
        #predicted_probability = open(fa_outname).readlines()[1]  
        predicted_secondary = parse_psipred_file(fa_filename+".ss")
        # remove files produced
        remove_file(fa_filename+".ss") 
        remove_file(fa_filename+".horiz") 
        remove_file(fa_filename+".ss2")
        # calculate metrics
        print(predicted_secondary)
        accuracy, confusion_m = prediction_rank(predicted_secondary,prot_true_sec)
        # update output dictionary
        output_dict[test_iteration][0].append(accuracy)
        output_dict[test_iteration][1].append(prot_test_seq_len)
        output_dict[test_iteration][2] += confusion_m

        # increment counter
        seq_count += 1
        
    return output_dict

def create_fasta(fa_filename, fasta_name, prot_test_seq):
    """
    This method creates a fasta file
    """
    file_out = open(fa_filename, "w+")
    file_out.write(f">{fasta_name} \n{prot_test_seq}")

    return None

def remove_file(file_to_delete):
    """
    This function removes a file on the system.
    """
    os.remove(file_to_delete)

    return None

def debruijnextend_test(python_path, train_pickle, test_df, test_iteration):
    """                                                                         
    this function performs the benchmarking for the software: psipredict        
                                                                                
    INPUT: 1. Path to the debruijnextend py file
           2. Path to the correct pickle file
           3. Test/true pandas dataframe
    OUTPUT:                                                                     
           dictionary = { 1 : [[.8, .4, ..., 0.7], # percent guessed correctly  
                               [254, 223, ..., 30], # length per                
                               [[], # confusion matrix                          
                                [],                                             
                                []                                              
                               ]                                                
                              ],                                                
                          ...                                                   
                        }                                                       
    """
    output_dict = {test_iteration : [[], [], np.zeros((3,3))] }                                                                     
    seq_count = 0
    for index, row in test_df.iterrows(): # loop through tests
        prot_test_seq = row['sequence']
        prot_test_seq_len = len(prot_test_seq)
        prot_true_sec = row['structure']
        # input for DebruijnExtend
        fa_filename = f"fa_filename_{test_iteration}_{seq_count}"
        fa_outname = f"fa_outname_{test_iteration}_{seq_count}"
        fasta_name = f"fasta_name_{test_iteration}_{seq_count}"
        # run the command
        create_fasta(fa_filename, fasta_name, prot_test_seq)
        CMD = f"python3 {python_path} {fa_filename} 4 {fa_outname}"
        print(CMD)
        os.system(CMD)
        remove_file(fa_filename)
        # parse the output
        predicted_probability = open(fa_outname).readlines()[1]  
        predicted_secondary = open(fa_outname).readlines()[2].strip('"').strip("'")
        remove_file(fa_outname)
        # calculate metrics
        accuracy, confusion_m = prediction_rank(predicted_secondary,prot_true_sec)
        # update output dictionary
        output_dict[test_iteration][0].append(accuracy)
        output_dict[test_iteration][1].append(prot_test_seq_len)
        output_dict[test_iteration][2] += confusion_m
        
        # increment counter
        seq_count += 1

    return output_dict

def save_as_pickle(data_structure, output_name):
    """
    This function is used to save the output benchmarking dictionaries as
    pickle files for ease of plotting. This way the benchmarking only has
    to happen once and plotting can happen seperately.
    """
    pickle.dump(data_structure, open( output_name, "wb" ) )


####
# MAIN FUNCTION
####    
def main():
    """
    This function controls the flow of the script as well as ensuring
    certain steps happen for benchmarking.
    """
    # Paths
    debruijn_path = "/Users/dreyceyalbin/Dropbox/Fall2020-classes/Algorithms/project/DebruijnExtend/"
    python_path = debruijn_path + "DebruijnExtend_v2.py"
    psipred_path = debruijn_path + "benchmarking/" + "runpsipred_single"
    kfold_dir = debruijn_path + "kfold_testandtrain/" 
    
    # init
    deb_dict = {}
    psipred_dict = {}

    # import the config file for testing and training
    benchmark_config = open(debruijn_path + "benchmarking/" + sys.argv[1]).readlines()
    
    # K-fold cross validation
    kfold_iteration = 0
    for line in benchmark_config:
        line = line.split(",")
        test = kfold_dir + line[0]
        train = kfold_dir + line[1]
        test_df = createpandas_df(test)
        # benchmark Debruijn-Extend
        deb_dict_temp = debruijnextend_test(python_path, 
                                            train,
                                            test_df,
                                            kfold_iteration)  
              
        deb_dict.update(deb_dict_temp)
        save_as_pickle(deb_dict, "deb_dict.p")
        # benchmark PsiPred
        psipred_dict_temp = benchmark_psipredict(psipred_path, 
                                                 test_df,
                                                 kfold_iteration)
        psipred_dict.update(psipred_dict_temp)
        save_as_pickle(psipred_dict, "psipred_dict.p")

        # increment counter
        kfold_iteration += 1

if __name__ == "__main__":
    main()
