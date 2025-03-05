#!/bin/bash

source config.env


echo "Checking R packages..."
Rscript -e "renv::status()"

tf=$(mktemp)

echo "Checking plink..."
if ./bin/plink2 --version > $tf 2>&1; then
    echo "All good!"
else
    cat $tf
fi

echo "Checking flashpca..."
if ./bin/flashpca --version > $tf 2>&1; then
    echo "All good!"
else
    cat $tf
fi

echo "Checking king..."
./bin/king > $tf 2>&1
if cat $tf | grep -q "KING"; then
    echo "All good!"
else
    cat $tf
fi

echo "Checking liftOver..."
gunzip -c bin/liftOver.gz > bin/liftOver
chmod 777 bin/liftOver
./bin/liftOver > $tf 2>&1
if cat $tf | grep -q "liftOver"; then
    echo "All good!"
else
    cat $tf
fi
rm bin/liftOver


echo "Checking cohort name..."
if [[ ! $cohort_name =~ ^[a-zA-Z0-9_]+$ ]]; then
    echo "Error: Please check your config.env. The variable cohort_name must be alphanumeric with no spaces. It will be used to generate your results files."
    exit 1
else
    echo "$cohort_name is a valid cohort name!"
fi


# check build = hg19 or hg38
if [[ $build != "hg19" && $build != "hg38" ]]; then
    echo "Error: build specified in config.env must be 'hg19' or 'hg38'"
    exit 1
fi


echo "Checking genotype input list..."
nchr=$(cat ${genotype_input_list} | grep -c '^')
# Check nchr = 22 or 23
if [[ $nchr -ne 22 && $nchr -ne 23 ]]; then
    echo "Error: Expected 22 or 23 chromosomes, but found $nchr"
    exit 1
fi

echo "Checking $nchr bgen files exist"
for i in $(seq 1 $nchr)
do
    bgen=$(awk -v i=$i 'NR==i { print $1 }' ${genotype_input_list})
    sample=$(awk -v i=$i 'NR==i { print $2 }' ${genotype_input_list})
    if [ ! -f "${bgen}" ]; then
        echo "${bgen} not found"
        exit 1
    fi

    if [ ! -f "${sample}" ]; then
        echo "${sample} not found"
        exit 1
    fi

    if [[ ! $bgen == *.bgen ]]
    then
        echo "$bgen should be a bgen file ending in .bgen"
        exit 1
    fi

    if [[ ! $sample == *.sample ]]
    then
        echo "$sample should be a sample file ending in .sample"
        exit 1
    fi
done
echo "All good!"
