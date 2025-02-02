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
md5sum -c ${cohort_name}_${step}.tar.gz.md5

sftp_address="lcgwassftp.blob.core.windows.net"
sftp_username_full="lcgwassftp.testcontainer.${sftp_username}"

read -s -p "Ready to upload? Press enter to continue: " anykey
echo ""
read -s -p "Enter SFTP password: " mypassword

export SSHPASS=$mypassword
${d}/bin/sshpass -e sftp -oBatchMode=no -b - ${sftp_username_full}@${sftp_address} << !
   put ${results_dir}/${cohort_name}_${step}.tar.gz
   put ${results_dir}/${cohort_name}_${step}.tar.gz.md5
   bye
!

# This used to work and now doesn't - need to figure it out
# curl -1 -v -k "sftp://${sftp_address}:22/${sftp_username}/${cohort_name}_${step}.tar.gz" --user "${sftp_username_full}:${mypassword}" -T "${results_dir}/${cohort_name}_${step}.tar.gz"
# curl -1 -v -k "sftp://${sftp_address}:22/${sftp_username}/${cohort_name}_${step}.tar.gz.md5" --user "${sftp_username_full}:${mypassword}" -T "${results_dir}/${cohort_name}_${step}.tar.gz.md5"

echo "Successfully uploaded ${cohort_name}_${step}.tar.gz and ${cohort_name}_${step}.tar.gz.md5. Many thanks indeed!"
