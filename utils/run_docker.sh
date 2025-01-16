#!/bin/bash

set -e

source config.env

# docker build -t mrcieu/lifecourse-gwas:latest .
# docker pull mrcieu/lifecourse-gwas:latest

genotype_input_dir=$(head -n 1 geno_input.txt | awk '{ print $1 }' | xargs -I {} dirname {})

docker run \
    -v "${PWD}/config.env:/project/config.env:ro" \
    -v "${PWD}/geno_input.txt:/project/geno_input.txt:ro" \
    -v ${phenotype_input_dir}:${phenotype_input_dir} \
    -v ${phenotype_processed_dir}:${phenotype_processed_dir} \
    -v ${genotype_input_dir}:${genotype_input_dir} \
    -v ${genotype_processed_dir}:${genotype_processed_dir} \
    -v ${results_dir}:${results_dir} \
    mrcieu/lifecourse-gwas:latest \
    "$@"
