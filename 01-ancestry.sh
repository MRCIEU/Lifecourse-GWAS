#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/01

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/01/logfile)

# Inputs:
# - Cleaned genotype data

# Processes:
# - Generate PCs (aware of relatedness if necessary)
# - Generate sparse GRM (for family data)

# Outputs:
# - Sparse GRM for each ancestry
# - PCs for each ancestry
# - mbfile for fastGWA

# Get bfile prefix names into an array
bfiles=( $(ls ${genotype_processed_dir}/bfiles/${bfile_prefix}*.bed | \
    xargs -I {} sh -c "basename {}" | \
    xargs -I {} sh -c "echo {} | sed 's/.bed//g'" ))

# bfiles=$(< ${genotype_processed_dir}/geno_chrs.txt)

for f in ${bfiles[@]}
do
    echo $f
done

# Create scratch directory
mkdir -p ${genotype_processed_dir}/scratch

# Get list of pruned SNPs
if test -f "resources/genotypes/${major_ancestry}_pruned_variants.txt.gz"; then
    echo "Found prune file"
    prunefile="${genotype_processed_dir}/scratch/indep.prune.in"
    gunzip -c resources/genotypes/${major_ancestry}_pruned_variants.txt.gz > ${prunefile}
fi

# Prune for PC etc
> ${genotype_processed_dir}/scratch/bfiles
for f in ${bfiles[@]}
do
    echo $f
    awk -f resources/genotypes/highldregionsb37.awk ${genotype_input_dir}/${f}.bim > ${genotype_processed_dir}/scratch/${f}_highldregions.txt

    if ! test -f "${genotype_processed_dir}/scratch/indep.prune.in"; then
        bin/plink2 \
            --threads ${env_threads} \
            --bfile ${genotype_input_dir}/${f} \
            --exclude ${genotype_processed_dir}/scratch/${f}_highldregions.txt \
            --indep-pairwise 10000 5 0.1 \
            --out ${genotype_processed_dir}/scratch/${f}_indep
        prunefile="${genotype_processed_dir}/scratch/${f}_indep.prune.in"
    fi

    bin/plink2 \
        --threads ${env_threads} \
        --bfile ${genotype_input_dir}/${f} \
        --extract ${prunefile} \
        --make-bed \
        --out ${genotype_processed_dir}/scratch/${f}_indep

    echo ${genotype_processed_dir}/scratch/${f}_indep >> ${genotype_processed_dir}/scratch/bfiles
done

# Generate single bfile
f1=`head -n 1 ${genotype_processed_dir}/scratch/bfiles`
sed -i 1d ${genotype_processed_dir}/scratch/bfiles
bin/plink2 \
    --threads ${env_threads} \
    --bfile $f1 \
    --pmerge-list bfile ${genotype_processed_dir}/scratch/bfiles \
    --make-bed \
    --out ${genotype_processed_dir}/scratch/indep \
    --maj-ref

# If family data
## get list of relateds and list of unrelateds
## generate pcs in unrelateds and project to relateds
## nobody is removed

# If not family data
## get list of relateds and list of unrelateds
## use list of unrelateds as keeplist going forwards


# Get relateds and unrelateds

bin/king \
    -b ${genotype_processed_dir}/scratch/indep.bed \
    --unrelated \
    --degree 3 \
    --cpus ${env_threads} \
    --prefix ${genotype_processed_dir}/scratch/king

cp ${genotype_processed_dir}/scratch/kingunrelated.txt ${genotype_processed_dir}/kingunrelated.txt

bin/plink2 \
    --threads ${env_threads} \
    --bfile ${genotype_processed_dir}/scratch/indep \
    --keep ${genotype_processed_dir}/scratch/kingunrelated.txt \
    --make-bed \
    --out ${genotype_processed_dir}/scratch/indep_unrelated

if [ "${env_family_data}" = "true" ]
then
    bin/plink2 \
        --threads ${env_threads} \
        --bfile ${genotype_processed_dir}/scratch/indep \
        --remove ${genotype_processed_dir}/scratch/kingunrelated.txt \
        --make-bed \
        --out ${genotype_processed_dir}/scratch/indep_related
fi

# Generate PCs

if test -f "${genotype_processed_dir}/${bfile_prefix}pc.txt"; then
    echo "pcafile already provided"
    Rscript resources/genotypes/genetic_outliers.r \
        ${genotype_processed_dir}/${bfile_prefix}pc.txt \
        ${env_pca_sd} \
        ${env_n_pcs} \
        ${genotype_processed_dir}/${bfile_prefix}_genetic_outliers.txt \
        ${results_dir}/01/pcaplot.png

    n_outliers=`wc -l ${genotype_processed_dir}/${bfile_prefix}_genetic_outliers.txt | awk '{ print $1 }'`
    if [ "${n_outliers}" = "0" ]; then
        echo "WARNING: there are $n_outliers genetic outliers based on the user-provided PCs"
        echo "We recommend one of the following"
        echo "- changing the env_pca_sd threshold"
        echo "- removing those outliers in ${genotype_processed_dir}/${bfile_prefix}_genetic_outliers.txt and recalculating the PCs"
        echo "- allowing the pipeline to calculate the PCs and remove outliers itself"
        exit 1
    fi
else
    if [ "${env_family_data}" = "true" ]
    then
        bin/king \
            -b ${genotype_processed_dir}/scratch/indep_unrelated.bed,${genotype_processed_dir}/scratch/indep_related.bed \
            --mds ${env_n_pcs} \
            --cpus ${env_threads} \
            --projection \
            --prefix ${genotype_processed_dir}/${bfile_prefix}
    else
        bin/king \
            -b ${genotype_processed_dir}/scratch/indep_unrelated.bed \
            --mds ${env_n_pcs} \
            --cpus ${env_threads} \
            --prefix ${genotype_processed_dir}/${bfile_prefix}
    fi

    # Remove unnecessary columns from pc file
    cut -d " " -f 1,2,7- ${genotype_processed_dir}/${bfile_prefix}pc.txt > ${genotype_processed_dir}/${bfile_prefix}pc.txt_formatted
    mv ${genotype_processed_dir}/${bfile_prefix}pc.txt_formatted ${genotype_processed_dir}/${bfile_prefix}pc.txt

    # Check PCs e.g. by plotting them
    Rscript resources/genotypes/genetic_outliers.r \
        ${genotype_processed_dir}/${bfile_prefix}pc.txt \
        ${env_pca_sd} \
        ${env_n_pcs} \
        ${genotype_processed_dir}/${bfile_prefix}_genetic_outliers.txt \
        ${results_dir}/01/pcaplot.png

    n_outliers=`wc -l ${genotype_processed_dir}/${bfile_prefix}_genetic_outliers.txt | awk '{ print $1 }'`

    if [ "${n_outliers}" = "0" ]
    then
        echo "No genetic outliers detected"
    else
        # Remove genetic outliers from data
        bin/plink2 \
            --threads ${env_threads} \
            --bfile ${genotype_processed_dir}/scratch/indep_unrelated \
            --remove ${genotype_processed_dir}/${bfile_prefix}_genetic_outliers.txt \
            --make-bed \
            --out ${genotype_processed_dir}/scratch/indep_unrelated

        if [ "${env_family_data}" = "true" ]
        then
            bin/plink2 \
                --threads ${env_threads} \
                --bfile ${genotype_processed_dir}/scratch/indep_related \
                --remove ${genotype_processed_dir}/${bfile_prefix}_genetic_outliers.txt \
                --make-bed \
                --out ${genotype_processed_dir}/scratch/indep_related

            bin/king \
                -b ${genotype_processed_dir}/scratch/indep_unrelated.bed,${genotype_processed_dir}/scratch/indep_related.bed \
                --mds ${env_n_pcs} \
                --cpus ${env_threads} \
                --projection \
                --prefix ${genotype_processed_dir}/${bfile_prefix}
        else
            bin/king \
                -b ${genotype_processed_dir}/scratch/indep_unrelated.bed \
                --mds ${env_n_pcs} \
                --cpus ${env_threads} \
                --prefix ${genotype_processed_dir}/${bfile_prefix}
        fi

        mv ${results_dir}/01/pcaplot.png ${results_dir}/01/pcaplot_round1.png


        # Remove unnecessary columns from pc file
        cut -d " " -f 1,2,7- ${genotype_processed_dir}/${bfile_prefix}pc.txt > ${genotype_processed_dir}/${bfile_prefix}pc.txt_formatted

        mv ${genotype_processed_dir}/${bfile_prefix}pc.txt_formatted ${genotype_processed_dir}/${bfile_prefix}pc.txt

        Rscript resources/genotypes/genetic_outliers.r \
            ${genotype_processed_dir}/${bfile_prefix}pc.txt \
            ${env_pca_sd} \
            ${env_n_pcs} \
            ${genotype_processed_dir}/${bfile_prefix}_genetic_outliers.txt \
            ${results_dir}/01/pcaplot.png
    fi
fi


# Generate sparse GRM

if [ "${env_family_data}" = "true" ]
then
    bin/king \
        -b ${genotype_processed_dir}/scratch/indep.bed \
        --related \
        --degree 3 \
        --cpus ${env_threads} \
        --prefix ${genotype_processed_dir}/scratch/king

    awk '{ print $1, $3, $14 }' ${genotype_processed_dir}/scratch/king.kin0 | grep -v "4th" | sed 1d > ${genotype_processed_dir}/scratch/king.kin0.formatted

    Rscript resources/genotypes/pedFAM.R \
        ${genotype_processed_dir}/scratch/indep.fam \
        ${genotype_processed_dir}/scratch/king.kin0.formatted \
        ${genotype_processed_dir}/${bfile_prefix}
fi

echo "Successfully generated PCs etc!"
