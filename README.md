![example workflow](https://github.com/Dreycey/DebruijnExtend/actions/workflows/github_actions.yml/badge.svg)

![Debruijn Extend](figures/debruijnextend_logo.png)

# DebruijnExtend

DebruijnExtend is a machine learning algorithm that constructs a probability-based de Bruijn graph using kmers of protein secondary structure fragments associated with a given primary sequence. Once created, edge contraction is used to find paths of highest probability to determine the overall predicted secondary structure.

# Installation

USING CONDA:

```
conda env create -f environment.yml
```

# Usage

Below are examples of general usage along with building the DB. The most detailed information on usage is featured on the [DebruijnExtend Wiki](https://github.com/Dreycey/DebruijnExtend/wiki).

## Quick Example (General usage)

-   run the following command to test dependencies (and general usage; kmer size of 4):

```
python DebruijnExtend.py --input examples/gfp.fasta --hash_table prot_hashtables/prothashtable_4.p --output_file example_out
```

## BUILD: Building the hashmap (custom kmer sizes)

-   To use a custom kmer size (here we use size of 6), run the following:

```
python3 DebruijnExtend.py --hash_table 6 -i examples/gfp.fasta -o result.ss3
```

This creates the following hash file: `HashtableData/hashtable_k6.pickle`

## RUN: different options for running the tool.

-   To run the tool on the custom kmer size, run the following:

```
python DebruijnExtend.py --input examples/gfp.fasta --hash_table HashtableData/hashtable_k6.pickle --output_file example_out
```

The input fasta file may be a single protein or a multifasta. The ordering of the output file will have the same format with the predicted 3-state secondary structure appearing below the name of the protein.
