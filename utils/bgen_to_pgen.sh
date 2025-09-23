#!/bin/bash

set -e

source config.env

pgendir=$1
inclusion_file=$2

if [ -z $pgendir ]; then
    echo "Usage: ./bgen_to_pgen.sh <pgendir>" <sample inclusion file>
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

    ./bin/plink2 --bgen ${bgen} ref-first --sample ${sample} --make-pgen --out ${pgendir}/${bn} --threads ${env_threads} --maf ${env_maf} --keep ${inclusion_file}
    echo "${pgendir}/${bn}" >> $tf
done

echo "Merging chromosomes"

./bin/plink2 --pfile ${pgendir}/chr1 --merge-list $tf --make-pgen --out ${pgendir}/allchr --threads ${env_threads}

echo "Removing per-chromosome pgen files"

for i in $(seq 1 ${nchr})
do
    rm ${pgendir}/${bn}.pgen
    rm ${pgendir}/${bn}.pvar
    rm ${pgendir}/${bn}.psam
done

echo "New pgen dataset: ${pgendir}/allchr.pgen, ${pgendir}/allchr.pvar, ${pgendir}/allchr.psam"
echo "Successfully converted to bgen files to pgen"
