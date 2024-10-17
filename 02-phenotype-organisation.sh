#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/02

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/02/logfile)

# Inputs:

# - Phenotype files
# - Covariates file
# - Ancestry ID lists
# - PCs per ancestry

# Processes:

# - Clean phenotypes
# - Generate summaries / distributions of phenotypes
# - Generate LD matrices for each phenotype subset
# - Generate PRS - phenotype associations for each phenotype subset

# Outputs:

# - Per time-point phenotype files
#   - /output/phenotypes/<phencode>_<agecode>_<ancestry>.txt
# - Covariate files
#   - /output/covariates/covs.txt - PCs 1-10 + sex
# - Summaries / distributions of phenotypes
#   - /results/02/<phen>_summary.rdata
# - LD matrices
# - PRS scores


echo "Organising phenotypes"
Rscript resources/render.r resources/phenotypes/organise_phenotypes.rmd ${results_dir}/02


echo "Generating list of phenotypes"
ls ${phenotype_processed_dir}/*.phen > ${phenotype_processed_dir}/phenolist
nphen=`cat ${phenotype_processed_dir}/phenolist | wc -l`
echo "Generated ${nphen} phenotype subsets"

echo "Generating LD matrices for each phenotype subset"

mk_phen_bfile () {
    ph=$1
    mkdir -p ${genotype_processed_dir}/scratch/tophits
    # Get rsids to keep
    awk '{ print $1 }' resources/genotypes/tophits/${ph}.txt > ${genotype_processed_dir}/scratch/tophits/${ph}.hits

    # Make tophits data
    > ${genotype_processed_dir}/scratch/tophits/${ph}.mergefile
    while read bf; do
        echo "${bf}"
        fn=$(basename -- ${bf})

        mkdir -p ${genotype_processed_dir}/scratch/tophits/temp
        n=$(grep -wf ${genotype_processed_dir}/scratch/tophits/${ph}.hits ${bf}.bim | wc -l)

        if [ "$n" -gt 0 ]; then
            bin/plink2 \
                --threads ${env_threads} \
                --bfile ${bf} \
                --extract ${genotype_processed_dir}/scratch/tophits/${ph}.hits \
                --exclude ${genotype_processed_dir}/bfiles/${fn}_vremove \
                --remove ${genotype_processed_dir}/bfiles/sremove \
                --make-bed \
                --out ${genotype_processed_dir}/scratch/tophits/temp/${fn}
            
            if test -f ${genotype_processed_dir}/scratch/tophits/temp/${fn}.bed; then
                echo "${genotype_processed_dir}/scratch/tophits/temp/${fn}" >> ${genotype_processed_dir}/scratch/tophits/${ph}.mergefile
            fi
        fi
    done < ${genotype_processed_dir}/geno_chrs.txt

    f1=`head -n 1  ${genotype_processed_dir}/scratch/tophits/${ph}.mergefile`
    sed -i 1d ${genotype_processed_dir}/scratch/tophits/${ph}.mergefile
    bin/plink2 \
        --threads ${env_threads} \
        --bfile $f1 \
        --pmerge-list bfile ${genotype_processed_dir}/scratch/tophits/${ph}.mergefile \
        --make-bed \
        --out ${genotype_processed_dir}/scratch/tophits/${ph}
    rm -r ${genotype_processed_dir}/scratch/tophits/temp

    # Generate score
    bin/plink2 \
        --bfile ${genotype_processed_dir}/scratch/tophits/${ph} \
        --score resources/genotypes/tophits/${ph}.txt \
        --out ${genotype_processed_dir}/scratch/tophits/${ph}
}

# phens=( $(cat ${phenotype_processed_dir}/phenolist | xargs -n1 basename | cut -d "_" -f 1 | sort | uniq) )
# echo "${phens[@]}"
# for phen in ${phens[@]}
# do
#     mk_phen_bfile $phen    
# done

# Correlation matrix for each phenotype
i=1
mkdir -p ${results_dir}/02/ldmat
phenolist=( $(cat ${phenotype_processed_dir}/phenolist) )
for phen in ${phenolist[@]}
do
    if [ -z $index ] || [ "$index" -eq "$i" ] ; then            
        filename=$(basename -- ${phen})
        filename="${filename%.*}"
        echo $filename
        # Get rsids to keep
        ph=$(echo $filename | cut -d "_" -f 1)
        echo $ph

        if [ ! -f ${genotype_processed_dir}/scratch/tophits/${ph}.bed ]; then
            mk_phen_bfile ${ph}
        fi

        mkdir -p ${genotype_processed_dir}/scratch/ldmats

        # Get IDs to keep
        awk '{ print $1, $2 }' ${phen} > ${genotype_processed_dir}/scratch/ldmats/keeptemp
        

        # Get LD matrix
        bin/plink2 \
            --threads ${env_threads} \
            --keep ${genotype_processed_dir}/scratch/ldmats/keeptemp \
            --bfile ${genotype_processed_dir}/scratch/tophits/${ph} \
            --r2-unphased ref-based bin4 yes-really \
            --out ${results_dir}/02/ldmat/${filename}

        rm -r ${genotype_processed_dir}/scratch/ldmats
    fi
    i=$((i+1))
done


echo "Generating PRS-phenotype associations for each subset"

Rscript resources/phenotypes/score.r ${phenotype_processed_dir}/phenolist ${genotype_processed_dir}/scratch/tophits ${results_dir}/02

echo "Successfully organised phenotypes!"
