#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/01

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/01/logfile)



echo "Get list of pruned SNPs"
if test -f "resources/genotypes/${major_ancestry}_pruned_variants.txt.gz"; then
    echo "Found prune file"
    prunefile="${genotype_processed_dir}/scratch/indep.prune.in"
    gunzip -c resources/genotypes/${major_ancestry}_pruned_variants.txt.gz > ${prunefile}
fi

# Get list of tophits

> ${genotype_processed_dir}/scratch/tophitsnps.txt
for f in resources/genotypes/tophits/*.txt
do
    awk '{ print $1 }' $f >> ${genotype_processed_dir}/scratch/tophitsnps.txt
done
cat ${genotype_processed_dir}/scratch/tophitsnps.txt | sort | uniq >> ${prunefile}
cut -d "_" -f 1 ${prunefile} | tr ":" " " | awk '{ print $1, $2, $2, $1":"$2 }' > ${prunefile}_range


bgen="/local-scratch/data/ukb/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr19.bgen"

mkdir -p ${genotype_processed_dir}/bgen_extract

> ${genotype_processed_dir}/bgen_extract/mergelist

for bgen in ${bgenfiles[@]}
do
    plink2 \
        --bgen ${bgen} ref-first \
        --sample ${samplefile} \
        --extract ${prunefile}_range \
        --make-bed \
        --out ${genotype_processed_dir}/bgen_extract/$(basename ${bgen} .bgen)
    echo "${genotype_processed_dir}/bgen_extract/$(basename ${bgen} .bgen)" >> ${genotype_processed_dir}/bgen_extract/mergelist
done

f1=`head -n 1 ${genotype_processed_dir}/bgen_extract/mergelist`
sed -i 1d ${genotype_processed_dir}/bgen_extract/mergelist
bin/plink2 \
    --threads ${env_threads} \
    --bfile $f1 \
    --pmerge-list bfile ${genotype_processed_dir}/bgen_extract/mergelist \
    --make-bed \
    --out ${genotype_processed_dir}/scratch/indep \
    --maj-ref

rm -r ${genotype_processed_dir}/bgen_extract

Rscript resources/genotypes/variant_ids.r ${genotype_processed_dir}/scratch/indep ${genotype_processed_dir}/scratch/indep bin/plink2 ${env_threads}

echo "Successfully extracted pruned variants from bgen files"
