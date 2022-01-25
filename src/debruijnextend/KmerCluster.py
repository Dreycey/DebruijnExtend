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
# in-house packages
from src.debruijnextend.utils import hamming_dist




@dataclass
class KmerCluster:
    clusters: Dict[str, List[str]]

    def get_close_kmers_clusters(self, hash_table, kmer, centroid_diff_threshold, top_N=1):
        """
        finds possible structures using clusers instead of all vs all
        """
        # find relasted clusters    
        kmers_to_look_at = []
        for centroid, cluster_kmers in tqdm(self.clusters.items()):
            if hamming_dist(centroid, kmer) < centroid_diff_threshold:
                kmers_to_look_at += [kmer_i for kmer_i in cluster_kmers]
        # use found kmers for further evaluation
        priority_queue = []
        highest_score = float("inf")
        for kmer_j in tqdm(kmers_to_look_at):
            hamming_score = hamming_dist(kmer_j, kmer)
            if hamming_score < highest_score:
                secondary_structs = hash_table[kmer_j]
                priority_queue.append((hamming_score, secondary_structs))
                priority_queue.sort(key=lambda a: a[0])
            if len(priority_queue) > top_N: priority_queue.pop(-1)
            highest_score = priority_queue[-1][0]
        # turn into output dictionary
        output_dict = {}
        for saved_res in priority_queue:
            output_dict.update(saved_res[1])
        return output_dict

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

