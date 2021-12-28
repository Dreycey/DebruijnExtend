# Phage Protein Analysis
This directory contains the code/info used to predict the secondary structures of phage proteins, and the code used to perform the analysis. 

## Steps

1. Download all Uniprot phage proteins
```
curl "https://www.uniprot.org/uniprot/?query=Phage&format=fasta&sort=score" > phage_proteins.fa;
```

2. Cluster the Uniprot phage proteins.
```
identityThresholds=(95 90 85 80 75 70 65)
for iThresh in ${identityThresholds[@]}; do
  echo "Working on Threshold $iThresh";
  cd-hit -i phage_proteins.fa -o phage_proteins_${iThresh}.fa -c ${iThresh} -T 3
done
```