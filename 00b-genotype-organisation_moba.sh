#!/bin/bash
 
# strict stop if there are any errors
set -e
 
# get environmental variables
source config.env
./bin/gcta64 \
  --bfile /home/grace.power/work/gpower/data/lifecourse_gwas_data_curation/geno/MoBaPsychGen_v1-ec-eur-batch-basic-qc \
  --pheno /home/grace.power/work/gpower/analysis/LifecourseGWAS_pipeline_run/genotype_input_dir/scratch/phenrand.txt \
  --fastGWA-lr \
  --extract /home/grace.power/work/gpower/analysis/LifecourseGWAS_pipeline_run/genotype_input_dir/variant_inclusion.txt \
  --keep /home/grace.power/work/gpower/analysis/LifecourseGWAS_pipeline_run/genotype_input_dir/sample_inclusion.txt \
  --thread-num 15 \
  --geno 0.1 \
  --maf 0.01 \
  --out /home/grace.power/work/gpower/analysis/LifecourseGWAS_pipeline_run/genotype_input_dir/scratch/phenrand

# Run null model summary script
Rscript resources/genotypes/nullgwas.r \
  ${genotype_processed_dir}/scratch/phenrand.fastGWA \
  ${results_dir}/00
 
echo "Successfully summarised and filtered variants"
 
