#!/bin/bash


##################################################################################################################################
###### 
###### Pipeline step:		2	
###### Start date		    17.10.2023
###### Update:              23.11.2023
###### Program:         	KING version 2.3.2 (https://www.kingrelatedness.com/Download.shtml)
###### Written by:		    Grace M. Power
###### Packages loaded:		N/A
###### Datafiles :		    [studyfile]_chr[chr].bed
######				        [studyfile]_chr[chr].fam
######				        [studyfile]_chr[chr].bim
######				        KGref.bed
######				        KGref.fam
######				        KGref.bim 
######				        (available at: https://www.kingrelatedness.com/ancestry/)
###### Objective :		    To generate labels for ancestrial group within each study sample using ancestry inference in KING. 
######                      One file is required as output to phenotype file development
######                      File per ancestry group are required to output to run further genotype analyses
######
##################################################################################################################################

# strict stop if there are any errors
set -e

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

source config.env

wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xzvf Linux-king.tar.gz

# Download the KING reference


mkdir -p ${genotype_input_dir}/king_test
cd king_test
mkdir out

wget -O ${genotype_input_dir}/king_test/KGref.bed.xz https://www.kingrelatedness.com/ancestry/KGref.bed.xz
wget -O ${genotype_input_dir}/king_test/KGref.fam.xz https://www.kingrelatedness.com/ancestry/KGref.fam.xz
wget -O ${genotype_input_dir}/king_test/KGref.bim.xz https://www.kingrelatedness.com/ancestry/KGref.bim.xz

# decompress files

xz -d -v ${genotype_input_dirK}/king_test/KGref.bed.xz 
xz -d -v ${genotype_input_dirK}/king_test/KGref.bim.xz
xz -d -v ${genotype_input_dirK}/king_test/KGref.fam.xz 

  ### update reference files to be reformated from rsid 

  ## chrom pos a1 a2

Rscript KGref

# 2. extract the king reference SNPs from our data (per chromosome)

awk '{print $2}' king_test/KGref.bim > king_test/KGref_snplist.txt

# 3. Run ancestry identification

module load languages/r/4.3.1
resource_dir=/user/work/sd20930/LifecourseMR_wg/ancestry_pipeline/king_test
source ${resource_dir}/config.sh
study_file=/mnt/storage/private/alspac/1/alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-06-25/data/genotypes/bestguess
study_name=ALSPAC
out_dir=/user/work/sd20930/LifecourseMR_wg/ancestry_pipeline/king_test/out

##

./king -b /user/work/sd20930/LifecourseMR_wg/ancestry_pipeline/king_test/KGref.bed,/mnt/storage/private/alspac/1/alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-06-25/data/genotypes/bestguess/data_chr21.bed --pca --degree --projection --rplot --prefix out_ancestry_ALSPAC

./king -b ${resource_dir}/KGref.bed, ${study_file}/data_chr${i}.bed --projection --pca --rplot --prefix ${out_dir}/out_ancestry_${study_name}


# TODO: Split the data into ancestry files, or have 1 file and have a ID keep list for each ancestry? Perhaps define a minimum sample size for an ancestral group
 

# 4. ld prune (per chromosome)
#variant pruning: window size in kb = 100, step size variant ct = 5, vif threshold = 1.05

for i in {01..22};
do
#plink --bfile ${genotype_input_dir}_chr${i} --extract ${resource_dir}/snplist.txt --prune parameters --make-bed --out ${genotype_input_dir}/${genotype_input_dir}${i}_extract

plink --bfile ${genotype_input_dir}/data_chr${i} -extract ${resource_dir}/KGref_snplist.txt -indep 100 5 1.05 -make-bed -out ${output_dir}/${i}_extract
done

### need to extract SNPs in {i}_extract.prune.in from {i}_extract.bed

# 5. combine into a single pruned .bed file

#plink --bmerge etc etc --make-bed --out ${genotype_input_dir}/${genotype_input_dir}_pruned
for i in {01..22};
do
plink --bfile data_chr${i} --exclude snps_to_exclude --make-bed --out data_chr${i}_pruned 

echo ${i}_extract >> mergelist.txt
done


plink -merge-list mergelist.txt -make-out -out merged

## Intermediate files: 
## OUTPUT: table with columns ID and ancestry

# TODO: check the syntax for jq
genotype_input_dir=$(jq "genotype_input_dir" config.json)
bfile_prefix=$(jq "bfile_prefix" config.json)

# TODO: Store logs somewhere that is being transferred back to HQ

# TODO: Allow people to just provide an ancestry mapping of IDs already, but specify the format
