#!/bin/bash
 
# strict stop if there are any errors
set -e
 
# get environmental variables
source config.env
./bin/gcta64 \
    --bed /home/grace.power/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.bed \
    --bim /home/grace.power/archive/moba/geno/MobaPsychgenReleaseMarch23/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.bim \
    --fam /home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/geno/MoBaPsychGen_v1-ec-eur-batch-basic-qc.fid_replaced.fam \   
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
 
