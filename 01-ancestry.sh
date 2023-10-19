#!/bin/bash

##################################################################################################################################
###### 
###### Pipeline step:		1	
###### Date			17.10.2023
###### Program:         	KING version 2.3.2 (https://www.kingrelatedness.com/Download.shtml)
###### Written by:		Grace M. Power
###### Packages loaded:		N/A
###### Datafiles :		[studyfile].bed
######				[studyfile].fam
######				[studyfile].bim
######				KGref.bed
######				KGref.fam
######				KGref.bim 
######				(available at: https://www.kingrelatedness.com/ancestry/)
###### Objective :		To generate ancestrial group labels within each study sample using ancestry inference in KING 
######
##################################################################################################################################

##[should I add a section of code to create plink binary files? These are required to run the below - [studyfile].bed etc] 
# https://zzz.bwh.harvard.edu/plink/data.shtml#bed

#In Linux download executable file and unzip

## TODO - define where we want to keep these downloads - probably in /resources
wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xzvf Linux-king.tar.gz


# 1. download the king reference
# 2. extract the king reference SNPs from our data (per chromosome)
# 3. ld prune (per chromosome)
# 4. combine into a single pruned file
# 5. Run ancestry identification


## Intermediate files: 
## OUTPUT: table with columns ID and ancestry

# get the input file names

# TODO: check the syntax for jq
genotype_input_dir=$(jq "genotype_input_dir" config.json)
bfile_prefix=$(jq "bfile_prefix" config.json)



# 2. extract the king reference SNPs from our data (per chromosome)

awk '{ print $2}' /path/to/king/king.bim > /path/to/king/snplist.txt

for i in {1..22};
do
    plink --bfile ${genotype_input_dir}/${genotype_input_dir}${i} --extract /path/to/king/snplist.txt --prune parameters --make-bed --out ${genotype_input_dir}/${genotype_input_dir}${i}_extract
done

plink --bfile ${genotype_input_dir}/${genotype_input_dir}X --extract /path/to/king/snplist.txt --prune parameters --make-bed --out ${genotype_input_dir}/${genotype_input_dir}X_extract

# merge the data into one
plink --bmerge etc etc --make-bed --out ${genotype_input_dir}/${genotype_input_dir}_pruned


/path/to/king -b KGref.bed, ${genotype_input_dir}/${genotype_input_dir}_pruned.bed --pca --projection --rplot --prefix out_ancestry_[studyname]

# TODO: extract king reference SNPs by chr:pos rather than rsid;
# or update reference to have chr:pos_a1_a2, in which case we extract by variant ID as normal

# TODO: the plink and king binaries will be stored in the repository somewhere, so when calling them use that path

# TODO: figure out the output of the ancestry inference - does it need to be edited for future use

# TODO: Split the data into ancestry files, or have 1 file and have a ID keep list for each ancestry? Perhaps define a minimum sample size for an ancestral group

# TODO: Store logs somewhere that is being transferred back to HQ

# TODO: Allow people to just provide an ancestry mapping of IDs already, but specify the format
