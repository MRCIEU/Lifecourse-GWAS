#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/04

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/04/logfile_aggregate)

nchr=$(cat ${genotype_input_list} | grep -c '^')
echo $nchr

echo ${results_dir}/04
echo ${phenotype_processed_dir}/phenolist

Rscript resources/genotypes/process_regenie_results.r ${results_dir}/04 $nchr ${phenotype_processed_dir}/phenolist

echo "Successfully aggregated results"
