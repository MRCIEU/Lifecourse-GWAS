#!/bin/bash

set -e

source config.env

newdir=$1
dosage_field=$2

if [ -z $newdir ]; then
    echo "Usage: ./convert_vcf_to_bgen.sh <newdir>"
    exit 1
fi

mkdir -p ${newdir}

nchr=$(cat ${genotype_input_list} | grep -c '^')
dn=$(head -n 1 ${genotype_input_list} | awk '{ print $1 }' | xargs dirname)

tf=$(mktemp)
for i in $(seq 1 ${nchr})
do
    vcf=$(awk -v i=$i 'NR==i { print $1 }' ${genotype_input_list})
    bn=$(basename $vcf)
    bn=$(echo $bn | sed 's/.vcf.gz//g')
    bn=$(echo $bn | sed 's/.vcf//g')
    echo "Processing ${vcf}"
    ./bin/plink2 --vcf ${vcf} vcf-dosage=$dosage_field --export bgen-1.3 id-paste=iid --out ${newdir}/${bn} --threads ${env_threads}
    ./bin/bgenix -g ${newdir}/${bn}.bgen -index -clobber
    echo "${newdir}/${bn}.bgen ${newdir}/${bn}.sample" >> $tf
done

cp ${genotype_input_list} ${genotype_input_list}.original
mv ${tf} ${genotype_input_list}

echo "Original bgen files are now listed in ${genotype_input_list}.original"
echo "New bgen files are now listed in ${genotype_input_list}"

echo "Successfully converted to bgen1.3 and indexed"
