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


phenolist=( $(cat ${phenotype_processed_dir}/phenolist) )

# Allow specific analysis to be run
# Can take any number between 1:ngwas where ngwas is the number of rows in ${phenotype_processed_dir}/phenolist
index=$1
nphen=`cat ${phenotype_processed_dir}/phenolist | wc -l`

if [ -z $index ]
then
    echo "Running all $nphen GWASs"
elif [ ! -z $index ]; then
    re='^[0-9]+$'
    if ! [[ $index =~ $re ]] ; then
        # check if $index is in the phenolist array
        if [[ " ${phenolist[@]} " =~ " ${index} " ]]; then
            echo "Running GWAS for phenotype $index"
        else
            echo "error: Index is not a number or a valid phenotype"
            echo "Usage: ${0} [index number]"
            exit 1
        fi
    else
        if [ "$index" -gt "$nphen" ] ; then
            echo "error: Index is larger than number of phenotypes"
            echo "Usage: ${0} [index number]"
            exit 1
        fi
        echo "Running $index of $nphen GWASs"
    fi
fi

echo $index

## TODO
# copy bim files over to results/04

# Do GWAS for each phenotype
i=1
for phen in ${phenolist[@]}
do
    if [ -z $index ] || [[ "$index" == "$phen" ]] || [[ "$index" == "$i" ]] ; then
        filename=$(basename -- ${phen})
        filename="${filename%.*}"
        echo $filename
        covs=$(echo $phen | sed 's/.phen$/.covs/1')
        echo $covs
        echo "0" > ${phen}.flag
        if [ "$env_family_data" == "true" ]
        then
            echo "family"
            ( ./bin/gcta-1.94.1  \
                --mbfile ${genotype_processed_dir}/geno_chrs.txt \
                --fastGWA-mlm \
                --grm-sparse ${genotype_processed_dir}/${bfile_prefix} \
                --exclude ${genotype_processed_dir}/bfiles/vremove \
                --keep ${genotype_processed_dir}/related_keep.txt \
                --pheno ${phen} \
                --qcovar ${covs} \
                --thread-num ${env_threads} \
                --maf 0 \
                --geno 1 \
                --out ${results_dir}/04/${filename} ) \
                || ( echo "1" > ${phen}.flag )
            flag=`cat ${phen}.flag`
            echo $flag
            if [ "$flag" -eq "1" ] ; then
                echo "LMM failed. Trying linear model using unrelateds only"
                ./bin/gcta-1.94.1 \
                    --mbfile ${genotype_processed_dir}/geno_chrs.txt \
                    --fastGWA-lr \
                    --exclude ${genotype_processed_dir}/bfiles/vremove \
                    --keep ${genotype_processed_dir}/unrelated_keep.txt \
                    --pheno ${phen} \
                    --qcovar ${covs} \
                    --thread-num ${env_threads} \
                    --maf 0 \
                    --geno 1 \
                    --out ${results_dir}/04/${filename}
            fi
        else
            echo "not family"
            ./bin/gcta-1.94.1 \
                --mbfile ${genotype_processed_dir}/geno_chrs.txt \
                --fastGWA-lr \
                --pheno ${phen} \
                --exclude ${genotype_processed_dir}/bfiles/vremove \
                --keep ${genotype_processed_dir}/related_keep.txt \
                --qcovar ${covs} \
                --thread-num ${env_threads} \
                --maf 0 \
                --geno 1 \
                --out ${results_dir}/04/${filename}
        fi
        # compress GWAS
        # keep only b, se because all other info is constant across GWASs
        echo "Compressing output..."
        Rscript resources/genotypes/compress_gwas.r ${results_dir}/04/${filename}.fastGWA ${results_dir}/00/variants.txt ${genotype_processed_dir}/bfiles/vremove
        rm ${results_dir}/04/${filename}.fastGWA
    fi
    i=$((i+1))
done

if [ -z $index ]
then
    echo "Successfully performed GWAS of all phenotypes!"
else
    echo "Successfully performed GWAS of phenotype $index!"
fi