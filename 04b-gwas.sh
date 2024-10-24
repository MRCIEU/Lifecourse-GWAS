#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/04

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/04/logfile_step2_${1})

# Inputs:

# - sparse GRM - 01
# - PCs - 01
# - genotype data - 00
# - phenotypes - 02
# - covariates - 00

# Processes:

# - FastGWA per phen x age x ancestry

# Output:

# - GWAS summary stats per phen x age x ancestry
#     - results/04/phen_<phencode>_<ancestry>_<age>.*


# Run step 1

chr=$1

# Check chr is between 1-23
if [[ $chr -lt 1 || $chr -gt 23 ]]; then
    echo "Error: Expected chromosome number between 1-23, but found $chr"
    echo "Usage: ./04b-gwas.sh <chr>"
    exit 1
fi


# Get list of bgen files
echo "Checking genotype input list..."
nchr=$(cat ${genotype_input_list} | grep -c '^')
# Check nchr = 22 or 23
if [[ $nchr -ne 22 && $nchr -ne 23 ]]; then
    echo "Error: Expected 22 or 23 chromosomes, but found $nchr"
    exit 1
fi

echo "Checking $nchr bgen files exist"
for i in $(seq 1 $nchr)
do
    bgen=$(awk -v i=$i 'NR==i { print $1 }' ${genotype_input_list})
    sample=$(awk -v i=$i 'NR==i { print $2 }' ${genotype_input_list})
    if [ ! -f "${bgen}" ]; then
        echo "${bgen} not found"
        exit 1
    fi

    if [ ! -f "${sample}" ]; then
        echo "${sample} not found"
        exit 1
    fi

    if [[ ! $bgen == *.bgen ]]
    then
        echo "$bgen should be a bgen file ending in .bgen"
        exit 1
    fi

    if [[ ! $sample == *.sample ]]
    then
        echo "$sample should be a sample file ending in .sample"
        exit 1
    fi
done
echo "All good!"


bgen=$(awk -v i=$chr 'NR==i { print $1 }' ${genotype_input_list})
sample=$(awk -v i=$chr 'NR==i { print $2 }' ${genotype_input_list})

# check if $sample is empty - this would mean it's a pgen fileset
if [ -z "$sample" ]; then
    bin/regenie_v3.6.gz_x86_64_Linux_mkl \
    --step 2 \
    --qt \
    --minINFO 0.8 \
    --pgen ${bgen} \
    --phenoFile ${phenotype_processed_dir}/regenie/phen.txt \
    --covarFile ${phenotype_processed_dir}/regenie/covs.txt \
    --bsize 1000 \
    --pred ${phenotype_processed_dir}/regenie/step1_pred.list \
    --out ${phenotype_processed_dir}/regenie/step2_${chr} \
    --threads ${env_threads}
else
    bin/regenie_v3.6.gz_x86_64_Linux_mkl \
    --step 2 \
    --qt \
    --minINFO 0.8 \
    --bgen ${bgen} \
    --sample ${sample} \
    --phenoFile ${phenotype_processed_dir}/regenie/phen.txt \
    --covarFile ${phenotype_processed_dir}/regenie/covs.txt \
    --bsize 1000 \
    --pred ${phenotype_processed_dir}/regenie/step1_pred.list \
    --out ${phenotype_processed_dir}/regenie/step2_${chr} \
    --threads ${env_threads}
fi

gzip ${phenotype_processed_dir}/regenie/step2_${chr}_*.regenie

echo "Successfully performed regenie step 2"
