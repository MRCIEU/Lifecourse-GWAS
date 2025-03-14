#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/00
mkdir -p ${genotype_processed_dir}/tmp
export TMPDIR=${genotype_processed_dir}/tmp


# log everything from this script to a logfile in the results directory
exec &> >(tee ${results_dir}/00/logfile_b)


echo "Get MAF and INFO scores"
Rscript resources/genotypes/rand.r ${genotype_processed_dir}/sample_inclusion.txt ${genotype_processed_dir}/scratch/phenrand.txt

samplefile=$(head -n 1 ${genotype_input_list} | awk '{print $2}')

./bin/gcta64 \
    --mbgen ${genotype_input_list} \
    --sample ${samplefile} \
    --pheno ${genotype_processed_dir}/scratch/phenrand.txt \
    --fastGWA-lr \
    --keep ${genotype_processed_dir}/sample_inclusion.txt \
    --thread-num ${env_threads} \
    --maf 0.005 \
    --geno 0.1 \
    --out ${genotype_processed_dir}/scratch/phenrand

Rscript resources/genotypes/organise_variants.r \
    ${genotype_processed_dir}/scratch/phenrand.fastGWA \
    ${results_dir}/00 \
    ${genotype_processed_dir}/variant_inclusion.txt \
    ${genotype_processed_dir}/build_mapping.txt


echo "Successfully summarised and filtered variants"
