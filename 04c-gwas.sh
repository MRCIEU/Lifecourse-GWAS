#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/04

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/04/logfile_aggregate)

nchr=$(cat ${genotype_input_list} | grep -c '^')
echo $nchr

Rscript ${results_dir}/04 $nchr ${phenotype_processed_dir}/phenolist



nphen=$(cat ${phenotype_processed_dir}/phenolist | grep -c '^')

cat ${phenotype_processed_dir}/phenolist | xargs basename

phenotype_processed_dir="/local-scratch/projects/Lifecourse-GWAS/gib/alspac/phen_proc2"
echo $phenotype_processed_dir

gwas=${phenotype_processed_dir}/$(cat ${phenotype_processed_dir}/phenolist | head -n 10 | tail -n 1)
echo $gwas
for gwas in $(cat ${phenotype_processed_dir}/phenolist)
do
    bn=$(basename $gwas | sed "s/.phen$//g")
    echo $bn
    out=${results_dir}/04/${bn}.regenie.gz
    > ${out}
    echo $out
    # for i in 1:nchr
    for i in $(seq 1 $nchr)
    do
        cat ${phenotype_processed_dir}/regenie/step2_${i}_${bn}.regenie.gz >> $out
    done    
done

ls -l /local-scratch/projects/Lifecourse-GWAS/gib/alspac/phen_proc2/regenie/step2_*_bmi_10-11_both.regenie.gz

ls -lh $out
ls -lh $out


less

libr


