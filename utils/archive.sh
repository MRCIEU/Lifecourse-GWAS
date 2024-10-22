#!/bin/bash

# strict stop if there are any errors
set -e

# get environmental variables
source config.env

step=$1
steps=(00 01 02 03 04)


if [ -z "$step" ]; then
  echo "Please provide a step number"
  echo "List of steps: ${steps[@]}"
  exit 1
fi



# Check if $step is in the list of steps
if [[ ! " ${steps[@]} " =~ " ${step} " ]]; then
  echo "Step $step is not in the list of steps"
  echo "List of steps: ${steps[@]}"
  exit 1
fi


echo "Checking that the step has completed"

if [ ! $step -eq "04" ]; then
    if $(tail -n 1 ${results_dir}/$step/logfile | grep -i -q "success") ;
    then
        echo "Step $step has completed successfully"
    else
        echo "Sorry, based on the logfile, step $step has not completed successfully. Please check the logfile or contact the developers for help."
        exit 1
    fi
fi

if [ $step -eq "04" ]; then
  
	nphen=`cat ${phenotype_processed_dir}/phenolist | wc -l`
	nsuccess=`tail -n 1 ${results_dir}/04/logfile* | grep -i "Success" | wc -l`
	if [ "${nphen}" = "${nsuccess}" ]; then
		echo "GWAS completed successfully for all ${nphen} phenotypes"
	else
		echo "Problem: only ${nsuccess} of ${nphen} expected phenotypes completed. Please check logs in ${results_dir}/04"
		exit 1
	fi
fi

echo "Tarballing results for step $step"

## Check that $cohort_name is alphanumeric with no spaces etc
if [[ ! $cohort_name =~ ^[a-zA-Z0-9_]+$ ]]; then
    echo "Error: Please check your config.env. The variable cohort_name must be alphanumeric with no spaces. It will be used to generate your results files."
    exit 1
fi

if [[ $step = "00" ]];
then
    cp config.env ${results_dir}/00
fi

tar -czvf ${results_dir}/${cohort_name}_$step.tar.gz -C ${results_dir} $step

d=$(pwd)
cd ${results_dir}
md5sum ${cohort_name}_$step.tar.gz > ${cohort_name}_$step.tar.gz.md5
cd $d

echo "Success: Results for step $step have been tarballed and checksummed. Thank you very much for running this step successfully!"
