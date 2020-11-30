#!/bin/bash

####
# This script is creates pickled dictionaries/hash tables for k-mers
# of sizes 3 to 31.
####


#######################################
# Make hash tables of multiple kmers sizes
# Arguments:
#   $1 CSV file with functions and structures
#   $2 The prefix for the output pickle file 
# Outputs:
#   This calls on the python object to create hash
#   tables of multiple sizes
#######################################
function make_hash_tables() {
    for k in {3..5}; do
        echo "python py_scripts/csvtohash.py $1 ${2}_${k}.p ${k}";
        python py_scripts/csvtohash.py $1 ${2}_${k}.p ${k}; 
    done
};

#######################################
# Main function controls the flow of the bash script.
# Arguments:
#   None
# Outputs:
#   STEP 1: makes hash tables of multiple sizes.
#######################################
function main() {
    # STEP 1: make the hash tables for each kmer size
    echo; echo "Now constructing hash tables.."; echo; echo;
    for csv_file in kfold_testandtrain/*train*;
    do
        make_hash_tables $csv_file ${csv_file%".csv"}_"prothash"
    done
}
main;
