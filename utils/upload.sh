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

d=$(pwd)
cd ${results_dir}
md5sum -c ${cohort_name}_${major_ancestry}_${step}.tar.gz.md5


${d}/bin/azcopy copy "${results_dir}/${cohort_name}_${major_ancestry}_${step}.tar.gz" "https://lcgwassftp.blob.core.windows.net/analysts/${sftp_username}/${cohort_name}_${major_ancestry}_${step}.tar.gz?${sas_token}"

${d}/bin/azcopy copy "${results_dir}/${cohort_name}_${major_ancestry}_${step}.tar.gz.md5" "https://lcgwassftp.blob.core.windows.net/analysts/${sftp_username}/${cohort_name}_${major_ancestry}_${step}.tar.gz.md5?${sas_token}"

echo "Successfully uploaded ${cohort_name}_${major_ancestry}_${step}.tar.gz and ${cohort_name}_${major_ancestry}_${step}.tar.gz.md5. Many thanks indeed!"
