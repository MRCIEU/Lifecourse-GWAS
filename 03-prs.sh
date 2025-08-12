#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/03
mkdir -p ${results_dir}/03/ldmats

mkdir -p ${genotype_processed_dir}/tmp
export TMPDIR=${genotype_processed_dir}/tmp

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/03/logfile)

echo "Updating plink file to have aligned effect alleles"
Rscript resources/genotypes/variant_ids.r ${genotype_processed_dir}/scratch/tophits ${genotype_processed_dir}/scratch/tophits2 bin/plink2 ${env_threads}

echo "Checking and removing duplicate variants from tophits2"
bin/plink2 \
  --bfile ${genotype_processed_dir}/scratch/tophits2 \
  --rm-dup force-first \
  --make-bed \
  --out ${genotype_processed_dir}/scratch/tophits2_nodups

mv ${genotype_processed_dir}/scratch/tophits2_nodups.bed ${genotype_processed_dir}/scratch/tophits2.bed
mv ${genotype_processed_dir}/scratch/tophits2_nodups.bim ${genotype_processed_dir}/scratch/tophits2.bim
mv ${genotype_processed_dir}/scratch/tophits2_nodups.fam ${genotype_processed_dir}/scratch/tophits2.fam

phenolist=( $(cat ${phenotype_processed_dir}/phenolist) )

# Allow specific analysis to be run
# Can take any number between 1:ngwas where ngwas is the number of rows in ${phenotype_processed_dir}/phenolist
index=$1
nphen=`cat ${phenotype_processed_dir}/phenolist | wc -l`

if [ -z $index ]
then
    echo "Running all $nphen GWASs"
elif [ ! -z $index ]; then
    re='^[0-9]+$'
    if ! [[ $index =~ $re ]] ; then
        # check if $index is in the phenolist array
        if [[ " ${phenolist[@]} " =~ " ${index} " ]]; then
            echo "Running GWAS for phenotype $index"
        else
            echo "error: Index is not a number or a valid phenotype"
            echo "Usage: ${0} [index number]"
            exit 1
        fi
    else
        if [ "$index" -gt "$nphen" ] ; then
            echo "error: Index is larger than number of phenotypes"
            echo "Usage: ${0} [index number]"
            exit 1
        fi
        echo "Running $index of $nphen GWASs"
    fi
fi

echo $index

ls ${phenotype_processed_dir}/*.phen > ${phenotype_processed_dir}/phenolist
nphen=`cat ${phenotype_processed_dir}/phenolist | wc -l`
echo "Generated ${nphen} phenotype subsets"

echo "Generating LD matrices for each phenotype subset"

# Correlation matrix for each phenotype
i=1
mkdir -p ${results_dir}/03/ldmats
mkdir -p ${genotype_processed_dir}/scratch/tophits
phenolist=( $(cat ${phenotype_processed_dir}/phenolist) )
for phen in ${phenolist[@]}
do
    if [ -z $index ] || [[ "$index" == "$phen" ]] || [[ "$index" == "$i" ]] ; then
        filename=$(basename -- ${phen})
        filename="${filename%.*}"
        echo $filename
        # Get rsids to keep
        ph=$(echo $filename | cut -d "_" -f 1)
        if [[ "$ph" == "bioavail" ]]; then
            ph="bioavail_testosterone"
        fi
        echo $ph

        if [ ! -f ${genotype_processed_dir}/scratch/tophits/${ph}.bed ]; then
            bin/plink2 \
                --bfile ${genotype_processed_dir}/scratch/tophits2 \
                --score resources/genotypes/tophits/${genome_build}/${ph}.txt \
                --threads ${env_threads} \
                --out ${genotype_processed_dir}/scratch/tophits/${ph}
        fi

        if [ ! -f ${genotype_processed_dir}/scratch/tophits/${ph}.hits ]; then
            awk '{ print $1 }' resources/genotypes/tophits/${genome_build}/${ph}.txt > ${genotype_processed_dir}/scratch/tophits/${ph}.hits
        fi

        mkdir -p ${genotype_processed_dir}/scratch/ldmats

        # Get IDs to keep
        awk '{ print $1, $2 }' ${phen} > ${genotype_processed_dir}/scratch/ldmats/keeptemp
        
        # Get LD matrix
        bin/plink2 \
            --threads ${env_threads} \
            --keep ${genotype_processed_dir}/scratch/ldmats/keeptemp \
            --bfile ${genotype_processed_dir}/scratch/tophits2 \
            --extract ${genotype_processed_dir}/scratch/tophits/${ph}.hits \
            --r-unphased ref-based bin4 yes-really \
            --out ${results_dir}/03/ldmats/${filename}

        rm -r ${genotype_processed_dir}/scratch/ldmats
    fi
    i=$((i+1))
done

echo "Successfully generated correlation matrices for each phenotype!"
echo "Generating PRS-phenotype associations for each subset"
Rscript resources/phenotypes/score.r ${phenotype_processed_dir}/phenolist ${genotype_processed_dir}/scratch/tophits ${results_dir}/03
Rscript resources/render.r resources/genotypes/prs.rmd ${results_dir}/03

echo "Successfully generated scores!"
