#! usr/bin/python3
# std pkgs
import sys
import pandas as pd
import random

"""
This script chunks the input CSV into different chunks for testing and training
depending upon the input k parameter.

USAGE:
   python kfoldprep.py <input CSV> <outfile pref> <k folds>
EXAMPLE:
python kfoldprep.py data/secondarystructure_COMP533.csv outfile_prefix 10
"""
                                                                            
def split_csv(csv_path, outprefix, k_folds):
    """                                                                     
    This method parses an input csv and send the information into           
    a seperate method which constructs the hash table.                      
    """                
    # initialize                                                     
    file_by_lines = open(csv_path).readlines()
    random.shuffle(file_by_lines)
    chunk_length = len(file_by_lines) // k_folds 
    chunk_i = 0
    k = 1

    # split into chunks
    print(f"file length: {len(file_by_lines)}\n")
    print(f"CHUNKS: \n")
    chunk_ranges = []
    for chunk in range(0, len(file_by_lines), chunk_length):
        if chunk_i == chunk:
            continue
        print(k, chunk_i, chunk)
        chunk_ranges.append([chunk_i, chunk]) 
        # increment
        chunk_i = chunk
        k += 1

    # use chunks to make new files
    iteration = 0
    for chunkrange_i in range(len(chunk_ranges)):
        test_file = open(outprefix+f"_test_{iteration}.csv","a+")
        train_file = open(outprefix+f"_train_{iteration}.csv","a+")
        for chunkrange_j in range(len(chunk_ranges)):
            start = chunk_ranges[chunkrange_j][0]
            end = chunk_ranges[chunkrange_j][1]
            if chunkrange_i == chunkrange_j:
                for line in file_by_lines[start:end]:
                    test_file.write(line)
            else:
                for line in file_by_lines[start:end]:
                    train_file.write(line)
        iteration += 1    

def main():
    """
    This script splits up the training CSV into train and test segments based
    on a specified k-fold number.
    """
    # init
    csv_path = sys.argv[1]
    outprefix = sys.argv[2]
    k_folds = int(sys.argv[3])
    # split up the CSV
    split_csv(csv_path, outprefix, k_folds)

if __name__ == "__main__":
    main()
