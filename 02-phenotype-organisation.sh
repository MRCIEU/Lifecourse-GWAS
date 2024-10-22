#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/02

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/02/logfile)

# Inputs:

# - Phenotype files
# - Covariates file
# - Ancestry ID lists
# - PCs per ancestry

# Processes:

# - Clean phenotypes
# - Generate summaries / distributions of phenotypes
# - Generate LD matrices for each phenotype subset
# - Generate PRS - phenotype associations for each phenotype subset

# Outputs:

# - Per time-point phenotype files
#   - /output/phenotypes/<phencode>_<agecode>_<ancestry>.txt
# - Covariate files
#   - /output/covariates/covs.txt - PCs 1-10 + sex
# - Summaries / distributions of phenotypes
#   - /results/02/<phen>_summary.rdata
# - html report

echo "Organising phenotypes"
Rscript resources/render.r resources/phenotypes/organise_phenotypes.rmd ${results_dir}/02

echo "Generating list of phenotypes"
ls ${phenotype_processed_dir}/*.phen > ${phenotype_processed_dir}/phenolist
nphen=`cat ${phenotype_processed_dir}/phenolist | wc -l`
echo "Generated ${nphen} phenotype subsets"

echo "Formatting for regenie"
Rscript resources/phenotypes/regenie.r ${phenotype_processed_dir}/phenolist

echo "Successfully organised phenotypes!"