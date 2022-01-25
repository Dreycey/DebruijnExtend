#! usr/bin/bash

#########
# This script runs all of the benchmarkings for a give input
# multi fasta file.
#########

function run_debext_woclusters() {
    echo; echo "+++++++++++++++"; echo "Running DebruijnExtend predictions"; echo "+++++++++++++++"; echo;
    local debext_path=$1 # ../../DebruijnExtend.py
    local tool_directory=$2
    local kmer_size=$3
    cd ${tool_directory}
    for fasta_file in ./*; do
        echo ${fasta_file};
        if [[ ! -f ${fasta_file%.fna}.ss3 ]]; then
            (time \
            python3 ${debext_path} -i ${fasta_file} \
                                   -ht ${kmer_size} \
                                   -o ${fasta_file%.fna}.ss3 \
            ) 2> ${fasta_file%.fna}_time.txt 
        fi
    done
    cd ../;
}

function run_debext() {
    echo; echo "+++++++++++++++"; echo "Running DebruijnExtend predictions"; echo "+++++++++++++++"; echo;
    local debext_path=$1 # ../../DebruijnExtend.py
    local tool_directory=$2
    local kmer_size=$3
    local cluster_file=$4
    cd ${tool_directory}
    for fasta_file in ./*; do
        echo ${fasta_file};
        if [[ ! -f ${fasta_file%.fna}.ss3 ]]; then
            (time \
            python3 ${debext_path} -i ${fasta_file} \
                                   -ht ${kmer_size} \
                                   -o ${fasta_file%.fna}.ss3 \
                                   --use_clusters ${cluster_file} \
            ) 2> ${fasta_file%.fna}_time.txt 
        fi
    done
    cd ../;
}

function make_fastarepo() {
    local multifasta=$1
    local new_directory=$2
    local suffix=$3
    echo; echo "Makingfasta repo: ${new_directory}"; echo
    make_directory ${new_directory};
    python3 multifasta2single.py ${multifasta} ${new_directory} ${suffix};
}

function make_directory() {
    local new_directory=$1
    echo; echo "Making directory: ${new_directory}"; echo
    if [[ ! -d ${new_directory} ]]; then
        mkdir ${new_directory};
    else
        rm -rf ${new_directory};
        mkdir ${new_directory};
    fi
}

function main() {
    local multifasta_testing=$1
    local directory_suffix=$2
    # make output diectories.
    # test out 1 by 1 hammng distance
    # for kmer_size in 6; do
    #     cluster_file_name='_hammingDistance'
    #     NEW_DIRECTORY=debext_"${directory_suffix}_${kmer_size}_${cluster_file_name}/"
    #     echo "Testing kmer size ${kmer_size} and cluster ${cluster_file_name}";
    #     make_fastarepo ${multifasta_testing} ${NEW_DIRECTORY} .fna
    #     run_debext_woclusters ../../../DebruijnExtend.py ${NEW_DIRECTORY} \
    #                           ../../../HashtableData/hashtable_k${kmer_size}.pickle;
    # done;
    # use multifasta and directories with each tool.
    for kmer_size in 6; do
        for cluster_file in ../../ClusterData/*k${kmer_size}*.pickle; do
            cluster_file_name=${cluster_file##*/}
            cluster_file_name=${cluster_file_name%.pickle}
            NEW_DIRECTORY=debext_"${directory_suffix}_${kmer_size}_${cluster_file_name}/"
            echo "Testing kmer size ${kmer_size} and cluster ${cluster_file_name}";
            make_fastarepo ${multifasta_testing} ${NEW_DIRECTORY} .fna
            run_debext ../../../DebruijnExtend.py ${NEW_DIRECTORY} \
                                                ../../../HashtableData/hashtable_k${kmer_size}.pickle \
                                                ../${cluster_file};
        done;
    done;
}
echo; echo "Running Benchmarking for DebruijnExtend"; echo; echo;
main $1 $2;
