#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

# create results directory
mkdir -p ${results_dir}/01

# log everything from this script to a logfile in the results director
exec &> >(tee ${results_dir}/01/logfile)

mkdir -p ${genotype_processed_dir}/tmp
export TMPDIR=${genotype_processed_dir}/tmp

# Inputs:
# - Cleaned genotype data

# Processes:
# - Generate PCs (aware of relatedness if necessary)
# - Generate sparse GRM (for family data)

# Outputs:
# - Sparse GRM for each ancestry
# - PCs for each ancestry
# - mbfile for fastGWA

containsElement () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo "There is no method for ${1}."
	echo "Please run:"
	echo "./01-check_data [arg]"
	echo "where arg is an optional argument that can be one of:"
	printf '%s\n' ${@:2}
	return 1
}

arg="all"
declare -a sections=('all' 'relateds' 'pcs' 'grm' 'keeplists')

if [ -n "${1}" ]; then
	arg="${1}"
	containsElement ${1} ${sections[@]}
fi

section_message () {

	echo "-----------------------------------------------"
	echo ""
	echo "$1 section"
	echo ""
	echo "to run this part on its own type:"
	echo "$ ./01-check_data.sh $1"
	echo ""
	echo "-----------------------------------------------"
	echo ""
	echo ""

}

# If family data
## get list of relateds and list of unrelateds
## generate pcs in unrelateds and project to relateds
## nobody is removed

# If not family data
## get list of relateds and list of unrelateds
## use list of unrelateds as keeplist going forwards


if [ "$arg" = "relateds" ] || [ "$arg" = "all" ]
then
    section_message "relateds"
    echo "Get relateds and unrelateds"

    bin/king \
        -b ${genotype_processed_dir}/scratch/indep.bed \
        --unrelated \
        --degree 3 \
        --cpus ${env_threads} \
        --prefix ${genotype_processed_dir}/scratch/king

    cp ${genotype_processed_dir}/scratch/kingunrelated.txt ${genotype_processed_dir}/kingunrelated.txt

    # Make sure tophits are removed
    bin/plink2 \
        --threads ${env_threads} \
        --bfile ${genotype_processed_dir}/scratch/indep \
        --exclude range ${thfile} \
        --keep ${genotype_processed_dir}/scratch/kingunrelated.txt \
        --make-bed \
        --out ${genotype_processed_dir}/scratch/indep_unrelated

    if [ "${env_family_data}" = "true" ]
    then
        bin/plink2 \
            --threads ${env_threads} \
            --exclude range ${thfile} \
            --bfile ${genotype_processed_dir}/scratch/indep \
            --remove ${genotype_processed_dir}/scratch/kingunrelated.txt \
            --make-bed \
            --out ${genotype_processed_dir}/scratch/indep_related
    fi

    # Create tophits file if it doesn't exist already
    if test -f "${genotype_processed_dir}/scratch/tophits.bed"; then
        echo "Generating tophits file"
        bin/plink2 \
            --threads ${env_threads} \
            --extract range ${thfile} \
            --bfile ${genotype_processed_dir}/scratch/indep \
            --make-bed \
            --out ${genotype_processed_dir}/scratch/tophits

fi

if [ "$arg" = "pcs" ] || [ "$arg" = "all" ]
then

    section_message "pcs"
    echo "Generate PCs"

    if test -f "${genotype_processed_dir}/pcs.txt"; then
        echo "pcafile already provided"
        Rscript resources/genotypes/genetic_outliers.r \
            ${genotype_processed_dir}/pcs.txt \
            ${env_pca_sd} \
            ${env_n_pcs} \
            ${genotype_processed_dir}/genetic_outliers.txt \
            ${results_dir}/01/pcaplot.png

        n_outliers=`wc -l ${genotype_processed_dir}/genetic_outliers.txt | awk '{ print $1 }'`
        if [ "${n_outliers}" != "0" ]; then
            echo "WARNING: there are $n_outliers genetic outliers based on the user-provided PCs"
            echo "We recommend one of the following"
            echo "- changing the env_pca_sd threshold"
            echo "- removing those outliers in ${genotype_processed_dir}/genetic_outliers.txt and recalculatingg the PCs"
            echo "- delete the ${genotype_processed_dir}/pcs.txt file and allow the pipeline to calculate the PCs and remove outliers itself"
            exit 1
        fi
        echo "Success - PCs already calculated and no outliers detected"

        if [ ! "$arg" = "all" ]; then
            exit 0
        fi
    fi

    pcs_unrelated () {
        bin/flashpca \
            --bfile ${genotype_processed_dir}/scratch/indep_unrelated \
            --ndim ${env_n_pcs} \
            --outpc ${genotype_processed_dir}/scratch/fastpca_pcs_unrelated.txt \
            --outload ${genotype_processed_dir}/scratch/fastpca_loadings.txt \
            --outmeans ${genotype_processed_dir}/scratch/fastpca_meansd.txt \
            --numthreads ${env_threads} \
            --outval ${genotype_processed_dir}/scratch/fastpca_eigenvalues.txt \
            --outvec ${genotype_processed_dir}/scratch/fastpca_eigenvectors.txt \
            --outpve ${genotype_processed_dir}/scratch/fastpca_pve.txt
    }

    pcs_related () {
        bin/flashpca \
            --bfile ${genotype_processed_dir}/scratch/indep_related \
            --project \
            --inmeansd ${genotype_processed_dir}/scratch/fastpca_meansd.txt \
            --outproj ${genotype_processed_dir}/scratch/fastpca_pcs_related.txt \
            --inload ${genotype_processed_dir}/scratch/fastpca_loadings.txt \
            --numthreads ${env_threads}
    }

    if [ "${env_family_data}" = "true" ]
    then
        pcs_unrelated
        pcs_related
        sed -i 1d ${genotype_processed_dir}/scratch/fastpca_pcs_related.txt
        cat ${genotype_processed_dir}/scratch/fastpca_pcs_unrelated.txt ${genotype_processed_dir}/scratch/fastpca_pcs_related.txt > ${genotype_processed_dir}/pcs.txt

    else
        pcs_unrelated
        cp ${genotype_processed_dir}/scratch/fastpca_pcs_unrelated.txt ${genotype_processed_dir}/pcs.txt
    fi

    echo "Check PCs e.g. by plotting them"
    Rscript resources/genotypes/genetic_outliers.r \
        ${genotype_processed_dir}/pcs.txt \
        ${env_pca_sd} \
        ${env_n_pcs} \
        ${genotype_processed_dir}/genetic_outliers.txt \
        ${results_dir}/01/pcaplot.png

    n_outliers=`wc -l ${genotype_processed_dir}/genetic_outliers.txt | awk '{ print $1 }'`

    if [ "${n_outliers}" = "0" ]
    then
        echo "No genetic outliers detected"
    else
        echo "Remove genetic outliers from data"
        echo "Found ${n_outliers}. Removing them and recalculating PCs"
        bin/plink2 \
            --threads ${env_threads} \
            --bfile ${genotype_processed_dir}/scratch/indep_unrelated \
            --remove ${genotype_processed_dir}/genetic_outliers.txt \
            --make-bed \
            --out ${genotype_processed_dir}/scratch/indep_unrelated

        if [ "${env_family_data}" = "true" ]
        then
            bin/plink2 \
                --threads ${env_threads} \
                --bfile ${genotype_processed_dir}/scratch/indep_related \
                --remove ${genotype_processed_dir}/genetic_outliers.txt \
                --make-bed \
                --out ${genotype_processed_dir}/scratch/indep_related

            pcs_unrelated
            pcs_related            
            sed -i 1d ${genotype_processed_dir}/scratch/fastpca_pcs_related.txt
            cat ${genotype_processed_dir}/scratch/fastpca_pcs_unrelated.txt ${genotype_processed_dir}/scratch/fastpca_pcs_related.txt > ${genotype_processed_dir}/pcs.txt

        else
            pcs_unrelated
            cp ${genotype_processed_dir}/scratch/fastpca_pcs_unrelated.txt ${genotype_processed_dir}/pcs.txt
        fi

        mv ${results_dir}/01/pcaplot.png ${results_dir}/01/pcaplot_round1.png

        Rscript resources/genotypes/genetic_outliers.r \
            ${genotype_processed_dir}/pcs.txt \
            ${env_pca_sd} \
            ${env_n_pcs} \
            ${genotype_processed_dir}/genetic_outliers.txt \
            ${results_dir}/01/pcaplot.png
    fi
fi

if [ "$arg" = "grm" ] || [ "$arg" = "all" ]
then
echo "Generate sparse GRM"

    section_message "grm"

    if [ "${env_family_data}" = "true" ]
    then
        bin/king \
            -b ${genotype_processed_dir}/scratch/indep.bed \
            --related \
            --degree 3 \
            --cpus ${env_threads} \
            --prefix ${genotype_processed_dir}/scratch/king

        awk '{ print $1, $3, $14 }' ${genotype_processed_dir}/scratch/king.kin0 | grep -v "4th" | sed 1d > ${genotype_processed_dir}/scratch/king.kin0.formatted

        Rscript resources/genotypes/pedFAM.R \
            ${genotype_processed_dir}/scratch/indep.fam \
            ${genotype_processed_dir}/scratch/king.kin0.formatted \
            ${genotype_processed_dir}/sparsegrm
    fi
fi

if [ "$arg" = "keeplists" ] || [ "$arg" = "all" ]
then
    section_message "keeplists"

    echo "Final keep lists"
    # Unrelateds
    # +kingunrelated.txt
    # -genetic_outliers.txt

    cat ${genotype_processed_dir}/kingunrelated.txt | \
        grep -vw -f ${genotype_processed_dir}/genetic_outliers.txt > \
        ${genotype_processed_dir}/unrelated_keep.txt

    nunrelated=$(cat ${genotype_processed_dir}/unrelated_keep.txt | wc -l)

    echo "N Unrelated: ${nunrelated}"

    # Relateds
    # +fam file (everyone)
    # -genetic_outliers.txt

    if [ "${env_family_data}" = "true" ]
    then
        awk '{ print $1"\t"$2 }' ${genotype_processed_dir}/scratch/indep.fam | \
        grep -vw -f ${genotype_processed_dir}/genetic_outliers.txt > \
        ${genotype_processed_dir}/related_keep.txt
    fi

    nrelated=$(cat ${genotype_processed_dir}/related_keep.txt | wc -l)

    echo "N Related: ${nrelated}"

fi

echo "Successfully generated PCs etc!"
