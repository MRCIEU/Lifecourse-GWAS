#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/04

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/04/logfile)


echo "Successfully generated correlation matrices for each phenotype!"
echo "Generating PRS-phenotype associations for each subset"
Rscript resources/phenotypes/score.r ${phenotype_processed_dir}/phenolist ${genotype_processed_dir}/scratch/tophits ${results_dir}/04
echo "Successfully generated scores!"
