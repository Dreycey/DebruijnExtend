#! usr/bin/bash

#########
# This script runs all of the benchmarkings for a give input
# multi fasta file.
#########

function run_jpred() {
    echo; echo "+++++++++++++++"; echo "Running JPRED predictions"; echo "+++++++++++++++"; echo;
    local tool_directory=$1
    local emailaddress=$2
    if [[ ! -f JPREDDONE.txt ]]; then
        python3 -m jpredapi submit --mode=batch \
                                   --format=fasta \
                                   --file=$tool_directory \
                                   --email=$emailaddress;
        touch JPREDDONE.txt;
    fi
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

function run_s4pred() {
    echo; echo "+++++++++++++++"; echo "Running S4PRED predictions"; echo "+++++++++++++++"; echo;
    local s4pred_runmodel=$1 # ./s4pred/run_model.py
    local tool_directory=$2
    cd ${tool_directory}
    for fasta_file in ./*; do
        echo ${fasta_file};
        if [[ ! -f ${fasta_file%.fna}.ss3 ]]; then
            (time \
                python ${s4pred_runmodel} --outfmt fas ${fasta_file} > ${fasta_file%.fna}.ss3
            ) 2> ${fasta_file%.fna}_time.txt 
        fi
    done
    cd ../;
}

function run_spider2() {
    echo; echo "+++++++++++++++"; echo "Running SPIDER2 predictions"; echo "+++++++++++++++"; echo;
    local run_local_path=$1 #SPIDER2/misc/run_local.sh
    local tool_directory=$2 # files with .seq fasta files
    cd ${tool_directory}
    for fasta_file in ./*seq; do
        echo ${fasta_file};
        echo bash $run_local_path ${fasta_file}; 
        if [[ ! -f ${fasta_file%.seq}.spd3 ]]; then
            (time \
                bash $run_local_path ${fasta_file};
            ) 2> ${fasta_file%.seq}_time.txt 
        fi
    done
    cd ../;
}

function run_porter5() {
    echo; echo "+++++++++++++++"; echo "Running PORTER5 predictions"; echo "+++++++++++++++"; echo;
    local porter5_path=$1 # Porter5/Porter5.py
    local tool_directory=$2 # files with .seq fasta files
    local threads=$3
    cd ${tool_directory}
    for fasta_file in ./*; do
        if [[ ! -f ${fasta_file%.fna}.ss3 ]]; then
            (time \
                python ${porter5_path} -i ${fasta_file} --fast --cpu ${threads};
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
    make_fastarepo ${multifasta_testing} s4pred_${directory_suffix}/ .seq
    make_fastarepo ${multifasta_testing} spider2_${directory_suffix}/ .seq
    make_fastarepo ${multifasta_testing} porter5_${directory_suffix}/ .fna
    make_fastarepo ${multifasta_testing} debext_${directory_suffix}/ .fna
    # use multifasta and directories with each tool.
    run_jpred ${multifasta_testing} albindreycey@gmail.com
    run_debext ../../../DebruijnExtend.py debext_${directory_suffix}/ ../TRAINING.pickle ../cluster_file.pickle 
    run_s4pred ../s4pred/run_model.py s4pred_${directory_suffix}/
    run_porter5 ../Porter5/Porter5.py porter5_${directory_suffix}/ 6
    run_spider2 ../SPIDER2/misc/run_local.sh spider2_${directory_suffix}/
}
echo; echo "Running Benchmarking for DebruijnExtend"; echo; echo;
main $1 $2;
