#!/usr/bin/env bash

#######################################
# This script uses CD-HIT to cluster an input 
# fasta file acccording on multiple identity 
# thresholds
#######################################
input_fasta=$1
identityThresholds=(0.95 0.90 0.85 0.80 0.75 0.70 0.65)
for iThresh in ${identityThresholds[@]}; do
  echo "Working on Threshold $iThresh";
  echo cd-hit -i ${input_fasta} -o phage_proteins_${iThresh}.fa -c ${iThresh} -T;
  cd-hit -i ${input_fasta} -o phage_proteins_${iThresh}.fa -c ${iThresh} -T 3 -M 4000;
done
