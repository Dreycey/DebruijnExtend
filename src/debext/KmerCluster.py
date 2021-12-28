"""
This module defines a datastructure for creating
kmers.
"""
from dataclasses import dataclass
import pickle
from typing import Dict, List, Optional, Union
from tqdm import tqdm
import os
from os.path import exists




def hamming_dist(k1, k2):
    val = 0
    for ind, char in enumerate(k1):
        if char != k2[ind]: val += 1
    return val

@dataclass
class KmerCluster:
    clusters: Dict[str, List[str]]

    @classmethod
    def init_struct(self, outputfile, cluster_file, threshold=6):
        """ create the data structure from a passed kmer dictionary """
        if exists(cluster_file): 
            with open(cluster_file, 'rb') as cls_file:
                return pickle.load(cls_file)

        clusters = {}
        hash_table = pickle.load(open(outputfile, "rb"))
        seqs = hash_table.keys()
        # greedy cluster
        counter = 0
        for prot_kmer in tqdm(seqs):
            cluster_found = False
            for centroid_kmer in clusters.keys():
                if hamming_dist(prot_kmer, centroid_kmer) <= threshold:
                    clusters[centroid_kmer].append(prot_kmer)
                    cluster_found = True
                    break
            if not cluster_found:
                clusters[prot_kmer] = [prot_kmer]
            counter += 1
            if (counter % 10000) == 0:
                print(f"number of centroids: {len(clusters.keys())}")
                print(f"number of total seqs: {counter}")
        print(len(clusters))

        # save the pickle
        with open(cluster_file, "wb") as outfile:
            pickle.dump(KmerCluster(clusters), outfile, protocol=pickle.HIGHEST_PROTOCOL)

        return KmerCluster(clusters)

