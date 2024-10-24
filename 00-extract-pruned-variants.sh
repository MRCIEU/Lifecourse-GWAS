#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/00

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/00/logfile)

mkdir -p ${genotype_processed_dir}/scratch

echo "Get list of pruned SNPs"
if test -f "resources/genotypes/hm3_prune_th_${build}.bed.gz"; then
    echo "Found prune file"
    prunefile="${genotype_processed_dir}/scratch/indep.prune.in"
    gunzip -c resources/genotypes/hm3_prune_th_${build}.bed.gz > ${prunefile}
else
    echo "Error: Prune file resources/genotypes/hm3_prune_th_${build}.bed.gz not found"
    exit 1
fi


# Get list of bgen files
echo "Checking genotype input list..."
nchr=$(cat ${genotype_input_list} | grep -c '^')
# Check nchr = 22 or 23
if [[ $nchr -ne 22 && $nchr -ne 23 ]]; then
    echo "Error: Expected 22 or 23 chromosomes, but found $nchr"
    exit 1
fi

echo "Checking $nchr bgen files exist"
for i in $(seq 1 $nchr)
do
    bgen=$(awk -v i=$i 'NR==i { print $1 }' ${genotype_input_list})
    sample=$(awk -v i=$i 'NR==i { print $2 }' ${genotype_input_list})
    if [ ! -f "${bgen}" ]; then
        echo "${bgen} not found"
        exit 1
    fi

    if [ ! -f "${sample}" ]; then
        echo "${sample} not found"
        exit 1
    fi

    if [[ ! $bgen == *.bgen ]]
    then
        echo "$bgen should be a bgen file ending in .bgen"
        exit 1
    fi

    if [[ ! $sample == *.sample ]]
    then
        echo "$sample should be a sample file ending in .sample"
        exit 1
    fi
done
echo "All good!"

mkdir -p ${genotype_processed_dir}/bgen_extract

> ${genotype_processed_dir}/bgen_extract/mergelist

for i in $(seq 1 $nchr)
do
    bgen=$(awk -v i=$i 'NR==i { print $1 }' ${genotype_input_list})
    sample=$(awk -v i=$i 'NR==i { print $2 }' ${genotype_input_list})
    # check if $sample is empty - this would mean it's a pgen fileset
    if [ -z "$sample" ]; then
        ./bin/plink2 \
            --bgen ${bgen} \
            --sample ${sample} \
            --extract range ${prunefile} \
            --make-bed \
            --out ${genotype_processed_dir}/bgen_extract/$(basename ${bgen} .bgen) \
            --threads ${env_threads}
        echo "${genotype_processed_dir}/bgen_extract/$(basename ${bgen})" >> ${genotype_processed_dir}/bgen_extract/mergelist
    else
        ./bin/plink2 \
            --bgen ${bgen} ref-first \
            --sample ${sample} \
            --extract range ${prunefile} \
            --make-bed \
            --out ${genotype_processed_dir}/bgen_extract/$(basename ${bgen} .bgen) \
            --threads ${env_threads}
        echo "${genotype_processed_dir}/bgen_extract/$(basename ${bgen} .bgen)" >> ${genotype_processed_dir}/bgen_extract/mergelist
done

./bin/plink2 \
    --threads ${env_threads} \
    --pmerge-list bfile ${genotype_processed_dir}/bgen_extract/mergelist \
    --make-bed \
    --max-alleles 2 \
    --out ${genotype_processed_dir}/scratch/indep \
    --maf 0.01 \
    --maj-ref

Rscript resources/genotypes/variant_ids_bim.r ${genotype_processed_dir}/scratch/indep

echo "Successfully extracted pruned variants from bgen files"
