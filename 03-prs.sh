#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/03

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/03/logfile)


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

ls ${phenotype_processed_dir}/*.phen > ${phenotype_processed_dir}/phenolist
nphen=`cat ${phenotype_processed_dir}/phenolist | wc -l`
echo "Generated ${nphen} phenotype subsets"

echo "Generating LD matrices for each phenotype subset"

# Correlation matrix for each phenotype
i=1
mkdir -p ${results_dir}/03
phenolist=( $(cat ${phenotype_processed_dir}/phenolist) )
for phen in ${phenolist[@]}
do
    if [ -z $index ] || [[ "$index" == "$phen" ]] || [[ "$index" == "$i" ]] ; then
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
            --out ${results_dir}/03/${filename}

        rm -r ${genotype_processed_dir}/scratch/ldmats
    fi
    i=$((i+1))
done

if [ -z $index ]
then
    echo "Successfully generated scores and ld matrices!"
else
    echo "Successfully generated score and ld matrix for $index!"
fi
