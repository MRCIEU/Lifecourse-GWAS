#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/03

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/03/logfile)


# Inputs:

# - sparse GRM
# - PCs
# - genotype data
# - phenotypes
# - covariates

# Processes:

# - FastGWA per phen x age x ancestry

# Output:

# - GWAS summary stats per phen x age x ancestry
#     - results/03/phen_<phencode>_<age>_<ancestry>.*

generate dummy files to practice
for anc in AFR EUR
do
    mkdir -p ${genotype_processed_dir}/$anc
    for i in {1..22}
    do
        touch ${genotype_processed_dir}/${anc}/${bfile_prefix}${i}.bed
        touch ${genotype_processed_dir}/${anc}/${bfile_prefix}${i}.bim
        touch ${genotype_processed_dir}/${anc}/${bfile_prefix}${i}.fam
    done
done

# Get the list of relevant ancestries
ancestries=( $(ls ${genotype_processed_dir}/*/*.bed | \
    xargs -I {} sh -c "dirname {}" | \
    xargs -I {} sh -c "basename {}" | \
    uniq | \
    xargs) )
echo $ancestries


# Make mbfile - list of all the per-chr bfiles for each ancestry
for anc in $ancestries
do
    echo $anc
    ls ${genotype_processed_dir}/${anc}/*.bed | sed 's/.bed//g' > ${genotype_processed_dir}/${anc}/geno_chrs.txt
    cat ${genotype_processed_dir}/${anc}/geno_chrs.txt
done


# Get list of phenotypes
ls ${phenotype_processed_dir}/phen_*.txt > ${phenotype_processed_dir}/phenolist

# Do GWAS for each phenotype
for anc in $ancestries
do
    phenolist=( $(grep $anc ${phenotype_processed_dir}/phenolist) )
    for phen in ${phenolist}
    do
        filename=$(basename -- ${phen})
        filename="${filename%.*}"
        echo $filename
        if [ "$env_family_data" == "true" ]
        then
            echo "family"
            ./bin/gcta-1.94.1 \
                --mbfile ${genotype_processed_dir}/${anc}/geno_chrs.txt \
                --fastGWA-mlm \
                --grm-sparse ${genotype_processed_dir}/${anc}/sp_grm \
                --pheno ${phenotype_input_dir}/ \
                --qcovar ${phenotype_input_dir}/covs_${anc}.txt \
                --thread-num ${env_threads} \
                --maf 0 \
                --geno 1 \
                --out ${results_dir}/03/${filename}
        elif
            echo "not family"
            ./bin/gcta-1.94.1 \
                --mbfile ${genotype_processed_dir}/${anc}/geno_chrs.txt \
                --fastGWA-lr \
                --pheno ${phenotype_input_dir}/ \
                --qcovar ${phenotype_input_dir}/covs_${anc}.txt \
                --thread-num ${env_threads} \
                --maf 0 \
                --geno 1 \
                --out ${results_dir}/03/${filename}
        fi
        # compress GWAS
        # keep only b, se, af, n because all other info is constant across GWASs
        # if space is a real issue could sacrifice af and n
        Rscript resources/genetics/compress_gwas.r ${filename}
    done
done

echo "Successfully performed GWAS of all phenotypes"



