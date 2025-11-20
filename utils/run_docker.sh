#!/bin/bash

set -e

source config.env

# docker build -t mrcieu/lifecourse-gwas:latest .
# docker pull mrcieu/lifecourse-gwas:latest

if [[ ! -z $pgen_input ]]; then
    echo "Using pgen input"
    # Check that $pgen_input is set
    pgen_input_dir=$(dirname $pgen_input)

    docker run \
        --env APPTAINER="true" \
        -v "${PWD}/config.env:/project/config.env:ro" \
        -v ${phenotype_input_dir}:${phenotype_input_dir} \
        -v ${phenotype_processed_dir}:${phenotype_processed_dir} \
        -v ${pgen_input_dir}:${pgen_input_dir} \
        -v ${genotype_processed_dir}:${genotype_processed_dir} \
        -v ${results_dir}:${results_dir} \
        mrcieu/lifecourse-gwas:latest \
        "$@"
else
    echo "Using bgen input"
    genotype_input_dir=$(head -n 1 geno_input.txt | awk '{ print $1 }' | xargs -I {} dirname {})

    docker run \
        --env APPTAINER="true" \
        -v "${PWD}/config.env:/project/config.env:ro" \
        -v "${PWD}/geno_input.txt:/project/geno_input.txt:ro" \
        -v ${phenotype_input_dir}:${phenotype_input_dir} \
        -v ${phenotype_processed_dir}:${phenotype_processed_dir} \
        -v ${genotype_input_dir}:${genotype_input_dir} \
        -v ${genotype_processed_dir}:${genotype_processed_dir} \
        -v ${results_dir}:${results_dir} \
        mrcieu/lifecourse-gwas:latest \
        "$@"
fi
