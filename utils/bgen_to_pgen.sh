#!/bin/bash

set -e

source config.env

pgendir=$1

if [ -z $pgendir ]; then
    echo "Usage: ./bgen_to_pgen.sh <pgendir>"
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

    ./bin/plink2 --bgen ${bgen} ref-first --sample ${sample} --make-pgen --out ${pgendir}/${bn} --threads ${env_threads}
    echo "${pgendir}/${bn}" >> $tf
done

cp ${genotype_input_list} ${genotype_input_list}.original
mv ${tf} ${genotype_input_list}

echo "Original bgen files are now listed in ${genotype_input_list}.original"
echo "New pgen files are now listed in ${genotype_input_list}"

echo "Successfully converted to bgen files to pgen"
