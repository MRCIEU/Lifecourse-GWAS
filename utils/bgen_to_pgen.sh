#!/bin/bash

set -e

source config.env

pgendir=$1
inclusion_file=$2

if [ -z $pgendir ]; then
    echo "Usage: ./bgen_to_pgen.sh <pgendir> <sample inclusion file>"
    exit 1
fi

mkdir -p ${pgendir}

nchr=$(cat ${genotype_input_list} | grep -c '^')
dn=$(head -n 1 ${genotype_input_list} | awk '{ print $1 }' | xargs dirname)
mkdir -p $dn/pgen

tf=$(mktemp)
for i in $(seq 1 ${nchr})
do
    bgen=$(awk -v i=$i 'NR==i { print $1 }' ${genotype_input_list})
    sample=$(awk -v i=$i 'NR==i { print $2 }' ${genotype_input_list})
    bn=$(basename $bgen .bgen)
    dn=$(dirname $bgen)

    ./bin/plink2 --bgen ${bgen} ref-first --sample ${sample} --make-pgen --out ${pgendir}/${bn} --threads ${env_threads} --maf ${env_minmaf} --keep ${inclusion_file} --geno ${env_miss} --mind ${env_imiss}
    echo "${pgendir}/${bn}" >> $tf
done

echo "Merging chromosomes"

./bin/plink2 --pmerge-list $tf --make-pgen --out ${pgendir}/allchr --threads ${env_threads}

echo "Removing 'chr' prefix from pvar file"

cp ${pgendir}/allchr.pvar ${pgendir}/allchr.pvar.original
sed -i 's/chr//g' ${pgendir}/allchr.pvar

echo "Chromosomes in original pvar file"
awk '{ print $1 }' ${pgendir}/allchr.pvar.original | sort | uniq -c 

echo "Chromosomes in new pvar file"
awk '{ print $1 }' ${pgendir}/allchr.pvar | sort | uniq -c 

echo "New pgen dataset: ${pgendir}/allchr.pgen, ${pgendir}/allchr.pvar, ${pgendir}/allchr.psam"
echo "Successfully converted to bgen files to pgen"
