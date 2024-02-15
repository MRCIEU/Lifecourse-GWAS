#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/00

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/00/logfile)


## Check that the software all runs



## Genotype cleaning

# - Require hg19
# - Rename all SNPs to be chr:pos_a1_a2
# - Remove SNPs with info < 0.5
# - Remove SNPs with MAF < 0.01
# - HW < 1e-7
# - split genotype data by chr if necessary


## Check phenotype data

# Get bfile prefix names into an array
bfiles=( $(ls ${genotype_input_dir}/${bfile_prefix}*.bed | \
    xargs -I {} sh -c "basename {}" | \
    xargs -I {} sh -c "echo {} | sed 's/.bed//g'" ))

mkdir -p ${genotype_processed_dir}/bfiles

echo "List of input bfiles:"
# Make mbfile - list of all the per-chr bfiles for each ancestry
> ${genotype_processed_dir}/geno_chrs.txt

for f in ${bfiles[@]}
do
    echo $f
    # Clean data
    bin/plink2 \
        --bfile ${genotype_input_dir}/${f} \
        --rm-dup force-first \
        --maf ${env_minmaf} \
        --hwe ${env_hwe} \
        --geno ${env_miss} \
        --mind ${snp_imiss} \
        --chr 1-23 \
        --make-bed \
        --out ${genotype_processed_dir}/bfiles/${f}_temp \
        --threads ${env_threads}
    
    # Update variant IDs and effect allele coding
    Rscript resources/genotypes/variant_ids.r ${genotype_processed_dir}/bfiles/${f}_temp ${genotype_processed_dir}/bfiles/${f} bin/plink2

    rm ${genotype_processed_dir}/bfiles/${f}_temp.bed
    rm ${genotype_processed_dir}/bfiles/${f}_temp.bim
    rm ${genotype_processed_dir}/bfiles/${f}_temp.fam

    echo "${genotype_processed_dir}/bfiles/${f}" >> ${genotype_processed_dir}/geno_chrs.txt
done
