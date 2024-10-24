#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/04

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/04/logfile_step1)

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

# For some reason regenie doesn't like variants with freq = 0.5
bin/plink2 --bfile ${genotype_processed_dir}/scratch/indep --freq --out ${genotype_processed_dir}/scratch/indep
awk '$6 == 0.5' ${genotype_processed_dir}/scratch/indep.afreq | awk '{print $2}' > ${genotype_processed_dir}/scratch/indep.remove


bin/regenie_v3.6.gz_x86_64_Linux_mkl \
  --step 1 \
  --bed ${genotype_processed_dir}/scratch/indep \
  --phenoFile ${phenotype_processed_dir}/regenie/phen.txt \
  --covarFile ${phenotype_processed_dir}/regenie/covs.txt \
  --bsize 1000 \
  --exclude ${genotype_processed_dir}/scratch/indep.remove \
  --lowmem \
  --lowmem-prefix ${phenotype_processed_dir}/regenie/tmp_rg \
  --out ${phenotype_processed_dir}/regenie/step1 \
  --threads ${env_threads} \
  --force-qt

echo "Successfully performed regenie step 1"
