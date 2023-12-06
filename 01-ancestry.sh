#!/bin/bash

##################################################################################################################################
###### 
###### Pipeline step:		2	
###### Start date		    17.10.2023
###### Update:              23.11.2023
###### Program:         	KING version 2.3.2 (https://www.kingrelatedness.com/Download.shtml)
###### Written by:		    Grace M. Power
###### Packages loaded:		N/A
###### Datafiles :		    [studyfile]_chr_[chr].bed
######				        [studyfile]_chr_[chr].fam
######				        [studyfile]_chr_[chr].bim
######				        KGref.bed
######				        KGref.fam
######				        KGref.bim 
######				        (available at: https://www.kingrelatedness.com/ancestry/)
###### Objective :		    To generate labels for ancestrial group within each study sample using ancestry inference in KING. 
######                      One file is required as output to phenotype file development
######                      File per ancestry group are required to output to run further genotype analyses
######
##################################################################################################################################


# Inputs:

# - Cleaned genotype data
# - King reference files

# Processes:

# - Identify ancestries (KING)
# - Keep set of ancestries with N > threshold 
# - Generate sparse GRM per ancestry 
# - Generate PCs per ancestry


# Outputs:

# - ID list for each ancestry
#   - /output/ancestries/<ancestry>_id.txt
# - Sparse GRM for each ancestry
#   - /output/grm/<ancestry>.*
# - PCs for each ancestry
#   - /output/pcs/<ancestry>.eigenvec


#Download executable file and unzip

wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xzvf Linux-king.tar.gz

# Download the KING reference

mkdir king_test
cd king_test
mkdir out

wget https://www.kingrelatedness.com/ancestry/KGref.bed.xz
wget https://www.kingrelatedness.com/ancestry/KGref.fam.xz
wget https://www.kingrelatedness.com/ancestry/KGref.bim.xz

# decompress files

xz -d -v KGref.bed.xz 
xz -d -v KGref.bim.xz 
xz -d -v KGref.fam.xz 

# 2. extract the king reference SNPs from our data (per chromosome)

awk '{print $2}' king_test/KGref.bim > king_test/KGref_snplist.txt

test_dir=/user/work/sd20930/LifecourseMR_wg/ancestry_pipeline/king_test
source ${test_dir}/config.sh

# 3. ld prune (per chromosome)
#variant pruning: window size in kb = 100, step size variant ct = 5, vif threshold = 1.05

for i in {01..22};
do
#plink --bfile ${genotype_input_dir}_chr_${i} --extract ${resource_dir}/snplist.txt --prune parameters --make-bed --out ${genotype_input_dir}/${genotype_input_dir}${i}_extract

plink --bfile ${genotype_input_dir}/data_chr${i} -extract ${resource_dir}/KGref_snplist.txt -indep 100 5 1.05 -make-bed -out ${output_dir}/${i}_extract
done

### need to extract SNPs in {i}_extract.prune.in from {i}_extract.bed

# 4. combine into a single pruned .bed file

#plink --bmerge etc etc --make-bed --out ${genotype_input_dir}/${genotype_input_dir}_pruned
for i in {01..22};
do
echo ${i}_extract >> mergelist.txt
done

### hello

plink -merge-list mergelist.txt -make-out -out merged

# 5. Run ancestry identification

./king -b ${resource_dir}/KGref.bed, ${output_dir}/merged.bed -projection -rplot -prefix out_ancestry_${study_name}

## Intermediate files: 
## OUTPUT: table with columns ID and ancestry

# get the input file names

# TODO: check the syntax for jq
genotype_input_dir=$(jq "genotype_input_dir" config.json)
bfile_prefix=$(jq "bfile_prefix" config.json)




# TODO: extract king reference SNPs by chr:pos rather than rsid;
# or update reference to have chr:pos_a1_a2, in which case we extract by variant ID as normal

# TODO: the plink and king binaries will be stored in the repository somewhere, so when calling them use that path

# TODO: figure out the output of the ancestry inference - does it need to be edited for future use

# TODO: Split the data into ancestry files, or have 1 file and have a ID keep list for each ancestry? Perhaps define a minimum sample size for an ancestral group

# TODO: Store logs somewhere that is being transferred back to HQ

# TODO: Allow people to just provide an ancestry mapping of IDs already, but specify the format
