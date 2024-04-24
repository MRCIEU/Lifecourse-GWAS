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
for f in ${bfiles[@]}
do
    echo $f
done

echo ""
echo "Cleaning each chromosome"

# Make mbfile - list of all the per-chr bfiles for each ancestry
# List of samples to remove
# List of variants to remove
> ${genotype_processed_dir}/geno_chrs.txt
> ${genotype_processed_dir}/bfiles/sremove
> ${genotype_processed_dir}/bfiles/vremove

for f in ${bfiles[@]}
do
    echo $f

    if test -f ${genotype_input_dir}/${f}.bim.orig; then
        echo "bim file previously updated. Restoring..."
        mv ${genotype_input_dir}/${f}.bim.orig ${genotype_input_dir}/${f}.bim
    fi

    # Update variant IDs and effect allele coding
    Rscript resources/genotypes/variant_ids_bim.r ${genotype_input_dir}/${f}
    echo "Updated variant IDs"

    cat ${genotype_input_dir}/${f}.bim | { grep "_duplicate" || true; } > ${genotype_processed_dir}/bfiles/${f}_temp_duplicate

    # Clean data
    bin/plink2 \
        --bfile ${genotype_input_dir}/${f} \
        --freq \
        --hardy \
        --missing \
        --out ${genotype_processed_dir}/bfiles/${f}_temp \
        --threads ${env_threads}

    awk -v maf=${env_minmaf} '($6 < maf || $6 > (1-maf)) {print $2}' ${genotype_processed_dir}/bfiles/${f}_temp.afreq > ${genotype_processed_dir}/bfiles/${f}_temp_mafsnps

    if test -f ${genotype_processed_dir}/bfiles/${f}_temp.hardy; then
        awk -v hwe=${env_hwe} '($10 < hwe) {print $2}' ${genotype_processed_dir}/bfiles/${f}_temp.hardy > ${genotype_processed_dir}/bfiles/${f}_temp_hardysnps
    fi

    if test -f ${genotype_processed_dir}/bfiles/${f}_temp.hardy.x; then
        awk -v hwe=${env_hwe} '($14 < hwe) {print $2}' ${genotype_processed_dir}/bfiles/${f}_temp.hardy.x > ${genotype_processed_dir}/bfiles/${f}_temp_hardysnps
    fi
    
    awk -v miss=${env_miss} '($5 > miss) {print $2}' ${genotype_processed_dir}/bfiles/${f}_temp.vmiss | { grep -wv "ID" || true; } > ${genotype_processed_dir}/bfiles/${f}_temp_vmiss
    awk -v miss=${env_imiss} '($5 > miss) {print $1, $2}' ${genotype_processed_dir}/bfiles/${f}_temp.smiss | { grep -wv "#FID IID" || true; } >> ${genotype_processed_dir}/bfiles/sremove

    cat ${genotype_processed_dir}/bfiles/${f}_temp_duplicate ${genotype_processed_dir}/bfiles/${f}_temp_mafsnps ${genotype_processed_dir}/bfiles/${f}_temp_hardysnps ${genotype_processed_dir}/bfiles/${f}_temp_vmiss | sort | uniq > ${genotype_processed_dir}/bfiles/${f}_vremove
    cat ${genotype_processed_dir}/bfiles/${f}_vremove >> ${genotype_processed_dir}/bfiles/vremove

    echo "${genotype_input_dir}/${f}" >> ${genotype_processed_dir}/geno_chrs.txt
    echo "Removing $(cat ${genotype_processed_dir}/bfiles/${f}_vremove | wc -l) variants"
done

sort ${genotype_processed_dir}/bfiles/sremove | uniq > temp
mv temp ${genotype_processed_dir}/bfiles/sremove
echo "Removing $(cat ${genotype_processed_dir}/bfiles/sremove | wc -l) individuals due to missing data"
echo "Removing $(cat ${genotype_processed_dir}/bfiles/vremove | wc -l) variants in total"


echo "Successfully completed 00-check-data.sh"
