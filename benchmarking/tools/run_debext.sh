#! usr/bin/bash

#########
# This script runs all of the benchmarkings for a give input
# multi fasta file.
# Example:
#     directory_suffix
#########

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
                                   #--use_clusters ${cluster_file} \
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
    echo ${identity_cutoff}
    # make output diectories.
    make_fastarepo ${multifasta_testing} debext_${directory_suffix}/ .fna
    #echo "../../../HashtableData_i${identity_cutoff}/hashtable_k6.pickle"
    # use multifasta and directories with each tool.
    run_debext ../../../DebruijnExtend.py debext_${directory_suffix}/ \
                                          ../../../HashtableData/hashtable_k6.pickle \
                                          #../../../ClusterData_i${identity_cutoff}/cluster_k6_CT4.pickle
}
echo; echo "Running Benchmarking for DebruijnExtend"; echo; echo;
main $1 $2;
