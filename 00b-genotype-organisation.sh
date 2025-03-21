#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/00
mkdir -p ${genotype_processed_dir}/tmp
export TMPDIR=${genotype_processed_dir}/tmp


# log everything from this script to a logfile in the results directory
exec &> >(tee ${results_dir}/00/logfile_b)


echo "Checking genotype input list..."
nchr=$(cat ${genotype_input_list} | grep -c '^')
# Check nchr = 22 or 23
if [[ $nchr -ne 22 && $nchr -ne 23 ]]; then
    echo "Error: Expected 22 or 23 chromosomes, but found $nchr"
    exit 1
fi

echo "Get MAF and INFO scores"

mkdir -p ${genotype_processed_dir}/info_scores

for i in $(seq 1 $nchr)
do
    bgen=$(awk -v i=$i 'NR==i { print $1 }' ${genotype_input_list})
    sample=$(awk -v i=$i 'NR==i { print $2 }' ${genotype_input_list})
    ./bin/plink2 \
        --bgen ${bgen} ref-first \
        --sample ${sample} \
        --keep ${genotype_processed_dir}/sample_inclusion.txt \
        --freq cols=chrom,pos,ref,alt,altfreq,nobs,machr2 \
        --maf ${env_minmaf} \
        --out ${genotype_processed_dir}/info_scores/$(basename ${bgen} .bgen) \
        --threads ${env_threads}
done

Rscript resources/genotypes/organise_variants.r \
    ${genotype_processed_dir}/info_scores \
    ${results_dir}/00 \
    ${genotype_processed_dir}/variant_inclusion.txt \
    ${genotype_processed_dir}/build_mapping.txt


echo "Null GWAS"

Rscript resources/genotypes/rand.r ${genotype_processed_dir}/sample_inclusion.txt ${genotype_processed_dir}/scratch/phenrand.txt

samplefile=$(head -n 1 ${genotype_input_list} | awk '{print $2}')

./bin/gcta64 \
    --mbgen ${genotype_input_list} \
    --sample ${samplefile} \
    --pheno ${genotype_processed_dir}/scratch/phenrand.txt \
    --fastGWA-lr \
    --extract ${genotype_processed_dir}/variant_inclusion.txt \
    --keep ${genotype_processed_dir}/sample_inclusion.txt \
    --thread-num ${env_threads} \
    --geno 0.1 \
    --maf ${env_minmaf} \
    --out ${genotype_processed_dir}/scratch/phenrand

Rscript resources/genotypes/nullgwas.r \
    ${genotype_processed_dir}/scratch/phenrand.fastGWA \
    ${results_dir}/00

echo "Successfully summarised and filtered variants"
