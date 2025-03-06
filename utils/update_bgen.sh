#!/bin/bash

set -e

source config.env

newdir=$1

if [ -z $newdir ]; then
    echo "Usage: ./update_bgen.sh <newdir>"
    exit 1
fi

mkdir -p ${newdir}

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
    ./bin/plink2 --bgen ${bgen} ref-first --sample ${sample} --export bgen-1.3 id-paste=iid --out ${newdir}/${bn} --threads ${env_threads}
    ./bin/bgenix -g ${newdir}/${bn}.bgen -index -clobber
    echo "${newdir}/${bn}.bgen ${newdir}/${bn}.sample" >> $tf
done

cp ${genotype_input_list} ${genotype_input_list}.original
mv ${tf} ${genotype_input_list}

echo "Original bgen files are now listed in ${genotype_input_list}.original"
echo "New bgen files are now listed in ${genotype_input_list}"

echo "Successfully converted to bgen1.3 and indexed"
