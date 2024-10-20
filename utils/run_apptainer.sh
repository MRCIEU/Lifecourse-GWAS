#!/bin/bash

source config.env

# apptainer pull -F docker://mrcieu/lifecourse-gwas:latest

# Run without pwd binding

# d=${PWD}
# echo $d
# td=$(mktemp -d)
# cd $td

apptainer run \
    --containall \
    --bind "${PWD}/config.env:/project/config.env" \
    --bind "${phenotype_input_dir}:${phenotype_input_dir}" \
    --bind "${phenotype_processed_dir}:${phenotype_processed_dir}" \
    --bind "${genotype_input_dir}:${genotype_input_dir}" \
    --bind "${genotype_processed_dir}:${genotype_processed_dir}" \
    --bind "${results_dir}:${results_dir}" \
    --cwd /project \
    ${PWD}/lifecourse-gwas_latest.sif \
    "$@"

# cd $d
