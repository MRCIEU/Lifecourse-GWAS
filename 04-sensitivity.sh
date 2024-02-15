#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/04

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/04/logfile${1})



# Allow specific analysis to be run
# Can take any number between 1:ngwas where ngwas is the number of rows in ${phenotype_processed_dir}/phenolist
index=$1
nphen=`cat ${phenotype_processed_dir}/phenolist | wc -l`


if [ -z $index ]
then
    echo "Running all $nphen GWASs"
else
    re='^[0-9]+$'
    if ! [[ $index =~ $re ]] ; then
        echo "error: Index variable is not a number"
        echo "Usage: ${0} [index number]"
        exit 1
    fi

    if [ "$index" -gt "$nphen" ] ; then
        echo "error: Index is larger than number of phenotypes"
        echo "Usage: ${0} [index number]"
        exit 1
    fi
    echo "Running $index of $nphen GWASs"
fi


# Correlation matrix for each phenotype
i=1
phenolist=( $(cat ${phenotype_processed_dir}/phenolist) )
for phen in ${phenolist[@]}
do
    if [ -z $index ] || [ "$index" -eq "$i" ] ; then            
        filename=$(basename -- ${phen})
        filename="${filename%.*}"
        echo $filename

        mkdir -p ${genotype_processed_dir}/scratch/ldmats

        # Get IDs to keep
        awk '{ print $1, $2 }' ${phen} > ${genotype_processed_dir}/scratch/ldmats/keeptemp

        # Get rsids to keep
        ph=$(echo $filename | cut -d "_" -f 1)
        echo $ph
        awk '{ print $1 }' resources/genotypes/tophits/${ph}.txt > ${genotype_processed_dir}/scratch/ldmats/hits.txt
        
        # Get LD matrix
        > ${genotype_processed_dir}/scratch/ldmats/mergefile
        while read bf; do
            echo "${bf}"
            fn=$(basename -- ${bf})
            bin/plink2 \
                --threads ${env_threads} \
                --keep ${genotype_processed_dir}/scratch/ldmats/keeptemp \
                --bfile ${bf} \
                --extract ${genotype_processed_dir}/scratch/ldmats/hits.txt \
                --make-bed \
                --out ${genotype_processed_dir}/scratch/ldmats/${fn}
            
            if test -f ${genotype_processed_dir}/scratch/ldmats/${fn}.bed; then
                echo "${genotype_processed_dir}/scratch/ldmats/${fn}" >> ${genotype_processed_dir}/scratch/ldmats/mergefile
            fi
        done < ${genotype_processed_dir}/geno_chrs.txt

        f1=`head -n 1  ${genotype_processed_dir}/scratch/ldmats/mergefile`
        sed -i 1d ${genotype_processed_dir}/scratch/ldmats/mergefile
        bin/plink2 \
            --threads ${env_threads} \
            --bfile $f1 \
            --pmerge-list bfile ${genotype_processed_dir}/scratch/ldmats/mergefile \
            --make-bed \
            --out ${genotype_processed_dir}/scratch/ldmats/mergefile

        bin/plink2 \
            --threads ${env_threads} \
            --bfile ${genotype_processed_dir}/scratch/ldmats/mergefile \
            --r2-phased bin4 yes-really \
            --out ${results_dir}/04/${filename}

        rm -r ${genotype_processed_dir}/scratch/ldmats
    fi
    i=$((i+1))
done

echo "Successfully performed GWAS of all phenotypes!"
