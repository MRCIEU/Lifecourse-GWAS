#!/bin/bash

# Inputs:

# - Phenotype files
# - Covariates file
# - Ancestry ID lists
# - PCs per ancestry

# Processes:

# - Clean phenotypes
# - Generate summaries / distributions of phenotypes

# Outputs:

# - Per time-point phenotype files
#   - /output/phenotypes/<phencode>_<agecode>_<ancestry>.txt
# - Covariate files
#   - /output/covariates/covs.txt - PCs 1-10 + sex
# - Summaries / distributions of phenotypes
#   - /results/02/<phen>_summary.rdata

# todo:
# check where the outputs all wanna go

Rscript resources/render.r resources/phenotypes/organise_phenotypes.rmd


