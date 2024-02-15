#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/03

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/03/logfile${1})


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
#     - results/03/phen_<phencode>_<ancestry>_<age>.*


# Allow specific analysis to be run
# Can take any number between 1:ngwas where ngwas is the number of rows in ${phenotype_processed_dir}/phenolist
index=$1
nphen=`cat ${phenotype_processed_dir}/phenolist | wc -l`


if [ -z $index ]
then
    echo "Running all $nphen GWASs"
else
    re='^[0-9]+$'
    if ! [[ $index =~ $re ]] ; then
        echo "error: Index variable is not a number"
        echo "Usage: ${0} [index number]"
        exit 1
    fi

    if [ "$index" -gt "$nphen" ] ; then
        echo "error: Index is larger than number of phenotypes"
        echo "Usage: ${0} [index number]"
        exit 1
    fi
    echo "Running $index of $nphen GWASs"
fi


## TODO
# copy bim files over to results/03

# Do GWAS for each phenotype
i=1
phenolist=( $(cat ${phenotype_processed_dir}/phenolist) )
for phen in ${phenolist[@]}
do
    if [ -z $index ] || [ "$index" -eq "$i" ] ; then            
        filename=$(basename -- ${phen})
        filename="${filename%.*}"
        echo $filename
        covs=$(echo $phen | sed 's/.phen$/.covs/1')
        echo $covs
        if [ "$env_family_data" == "true" ]
        then
            echo "family"
            ./bin/gcta-1.94.1 \
                --mbfile ${genotype_processed_dir}/geno_chrs.txt \
                --fastGWA-mlm \
                --grm-sparse ${genotype_processed_dir}/${bfile_prefix} \
                --pheno ${phen} \
                --qcovar ${covs} \
                --thread-num ${env_threads} \
                --maf 0 \
                --geno 1 \
                --out ${results_dir}/03/${filename}
        else
            echo "not family"
            ./bin/gcta-1.94.1 \
                --mbfile ${genotype_processed_dir}/geno_chrs.txt \
                --fastGWA-lr \
                --pheno ${phen} \
                --qcovar ${covs} \
                --thread-num ${env_threads} \
                --maf 0 \
                --geno 1 \
                --out ${results_dir}/03/${filename}
        fi
        # compress GWAS
        # keep only b, se, af, n because all other info is constant across GWASs
        # if space is a real issue could sacrifice af and n
        # add variantid
        Rscript resources/genotypes/compress_gwas.r ${results_dir}/03/${filename}.fastGWA
        rm ${results_dir}/03/${filename}.fastGWA
    fi
    i=$((i+1))
done

echo "Successfully performed GWAS of all phenotypes!"
