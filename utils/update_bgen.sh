#!/bin/bash

set -e

source config.env

nchr=$(cat ${genotype_input_list} | grep -c '^')
dn=$(head -n 1 ${genotype_input_list} | awk '{ print $1 }' | xargs dirname)
mkdir -p $dn/bgen1.2

tf=$(mktemp)
for i in $(seq 1 ${nchr})
do
    bgen=$(awk -v i=$i 'NR==i { print $1 }' ${genotype_input_list})
    sample=$(awk -v i=$i 'NR==i { print $2 }' ${genotype_input_list})
    bn=$(basename $bgen .bgen)
    dn=$(dirname $bgen)

    ./bin/plink2 --bgen ${bgen} ref-first --sample ${sample} --export bgen-1.2 --out ${dn}/bgen1.2/${bn}
    ./bin/bgenix -g ${dn}/bgen1.2/${bn}.bgen -index -clobber
    echo "${dn}/bgen1.2/${bn}.bgen ${dn}/bgen1.2/${bn}.sample" >> $tf
done

cp ${genotype_input_list} ${genotype_input_list}.original
mv ${tf} ${genotype_input_list}

echo "Original bgen files are now listed in ${genotype_input_list}.original"
echo "New bgen files are now listed in ${genotype_input_list}"

echo "Successfully converted to bgen1.2 and indexed"