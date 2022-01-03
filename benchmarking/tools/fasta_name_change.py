#! usr/bin/python3
"""
This script changes a fasta name to 
an ordered number sequence for use in benchmarking.

REASON: this is needed for the benchmarking since
the file names are saved as the sequence name, which 
may include spaces in certain scenarios.
"""
import sys

if __name__ == "__main__":
    protein_id = 0
    with open(sys.argv[1], "r") as fasta_file:
        for line in fasta_file.readlines():
            line = line.strip("\n")
            if line[0] == ">":
                print(f">{protein_id}")
                protein_id += 1
            else:
                print(line)
