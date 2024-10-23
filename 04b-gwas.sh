#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/04

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/04/logfile${1})

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

bfiles=( $(ls ${genotype_processed_dir}/symlinks/${bfile_prefix}*.bed | \
    xargs -I {} sh -c "basename {}" | \
    xargs -I {} sh -c "echo {} | sed 's/.bed//g'" ))

bfile=${bfiles[$((${chr}-1))]}
echo $bfile


# Calculate 0.01 * 50000
nid=$(cat ${genotype_processed_dir}/symlinks/${bfile}.fam | wc -l)
minMAC=$(($nid/100))
echo $minMAC

bin/regenie_v3.6.gz_x86_64_Linux_mkl \
  --step 2 \
  --qt \
  --minINFO 0.8 \
  --minMAC $minMAC \
  --bed ${genotype_processed_dir}/symlinks/${bfile} \
  --phenoFile ${phenotype_processed_dir}/regenie_test/phen.txt \
  --covarFile ${phenotype_processed_dir}/regenie_test/covs.txt \
  --exclude ${genotype_processed_dir}/bfiles/vremove \
  --remove ${genotype_processed_dir}/bfiles/sremove \
  --bsize 200 \
  --pred ${phenotype_processed_dir}/regenie/step1_pred.list \
  --out ${phenotype_processed_dir}/regenie/step2_${chr} \
  --threads ${env_threads}

gzip ${phenotype_processed_dir}/regenie/step2_${chr}_*.regenie

echo "Successfully performed regenie step 2"
