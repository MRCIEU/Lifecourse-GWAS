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

# Cap threads for memory safety unless user explicitly overrides with KING_CPUS
KING_CPUS="${KING_CPUS:-${env_threads:-1}}"
if [[ "${KING_CPUS}" -gt 4 ]]; then
  echo "Capping KING threads from ${KING_CPUS} to 4 to reduce RAM pressure."
  KING_CPUS=4
fi

# Ensure inputs exist early
for ext in bed bim fam; do
  if [[ ! -s "${genotype_processed_dir}/scratch/indep.${ext}" ]]; then
    echo "ERROR: Missing ${genotype_processed_dir}/scratch/indep.${ext}"
    exit 2
  fi
done

# Helper: run KING unrelated with graceful fallback to plink2 --king-cutoff on OOM
run_king_unrelated() {
  local bed="${genotype_processed_dir}/scratch/indep.bed"
  local prefix="${genotype_processed_dir}/scratch/king"
  echo "Running KING (unrelated, degree 3) with ${KING_CPUS} threads..."
  set +e
  /usr/bin/time -v bin/king -b "${bed}" --unrelated --degree 3 --cpus "${KING_CPUS}" --prefix "${prefix}"
  rc=$?
  set -e
  if [[ $rc -eq 0 && -s "${prefix}unrelated.txt" ]]; then
    echo "KING unrelated completed."
    return 0
  fi
  echo "WARNING: KING failed (rc=${rc}), falling back to plink2 --king-cutoff"
  bin/plink2 \
    --threads "${env_threads:-1}" \
    --bfile "${genotype_processed_dir}/scratch/indep" \
    --king-cutoff 0.0884 \
    --out "${genotype_processed_dir}/scratch/p2_king"
  cp "${genotype_processed_dir}/scratch/p2_king.king.cutoff.in.id" "${prefix}unrelated.txt"
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

    # Use the robust helper (KING with fallback to plink2 --king-cutoff)
    run_king_unrelated

    cp ${genotype_processed_dir}/scratch/kingunrelated.txt ${genotype_processed_dir}/kingunrelated.txt

    # Make sure tophits are removed
    thfile="${genotype_processed_dir}/scratch/th.txt"
    gunzip -c resources/genotypes/tophitsnps_${genome_build}.bed.gz > ${thfile}

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
    if ! [ -e "${genotype_processed_dir}/scratch/tophits.bed" ]; then
        echo "Generating tophits file"
        bin/plink2 \
            --threads ${env_threads} \
            --extract range ${thfile} \
            --bfile ${genotype_processed_dir}/scratch/indep \
            --make-bed \
            --out ${genotype_processed_dir}/scratch/tophits
    fi
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
    else
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



        # Make sure tophits are removed
        # Create tophits file if it doesn't exist already
        if ! [ -e "${genotype_processed_dir}/scratch/tophits.bed" ]; then

            thfile="${genotype_processed_dir}/scratch/th.txt"
            gunzip -c resources/genotypes/tophitsnps_${genome_build}.bed.gz > ${thfile}

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

            echo "Generating tophits file"
            bin/plink2 \
                --threads ${env_threads} \
                --extract range ${thfile} \
                --bfile ${genotype_processed_dir}/scratch/indep \
                --make-bed \
                --out ${genotype_processed_dir}/scratch/tophits
        fi

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
            --cpus ${KING_CPUS} \
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

    # Make sure genetic_outliers.txt exists (empty if none)
    : > "${genotype_processed_dir}/genetic_outliers.txt"

    # Unrelateds = kingunrelated minus outliers (if any)
    if [ -s "${genotype_processed_dir}/genetic_outliers.txt" ]; then
        grep -vw -f "${genotype_processed_dir}/genetic_outliers.txt" \
            "${genotype_processed_dir}/kingunrelated.txt" \
            > "${genotype_processed_dir}/unrelated_keep.txt"
    else
        cp "${genotype_processed_dir}/kingunrelated.txt" \
           "${genotype_processed_dir}/unrelated_keep.txt"
    fi

    nunrelated=$(wc -l < "${genotype_processed_dir}/unrelated_keep.txt")
    echo "N Unrelated: ${nunrelated}"

    # Relateds only if family data
    if [ "${env_family_data}" = "true" ]; then
        awk '{ print $1"\t"$2 }' "${genotype_processed_dir}/scratch/indep.fam" | \
        ( [ -s "${genotype_processed_dir}/genetic_outliers.txt" ] \
          && grep -vw -f "${genotype_processed_dir}/genetic_outliers.txt" || cat ) \
        > "${genotype_processed_dir}/related_keep.txt"

        nrelated=$(wc -l < "${genotype_processed_dir}/related_keep.txt")
        echo "N Related: ${nrelated}"
    else
        echo "Family data = false; skipping related_keep.txt"
    fi
fi


echo "Successfully generated PCs etc!"
