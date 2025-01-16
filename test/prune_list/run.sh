#!/bin/bash


# Get 1kg in plink format
wget http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz
tar xzvf data_maf0.01_rs_ref.tgz

# get hm3 list
wget https://www.broadinstitute.org/files/shared/mpg/hapmap3/hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2
bzip2 -d hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map.bz2

awk '{ print $2 }' hapmap3_r1_b36_fwd_consensus.qc.poly.recode.map > hapmap3_r1_b36_fwd_consensus.qc.poly.recode.snp

# prune hm3 using 1kg data
plink2 --bfile data_maf0.01_rs_ref --exclude range excl_hg19.txt --out indep --make-bed

plink2 --bfile indep --extract hapmap3_r1_b36_fwd_consensus.qc.poly.recode.snp --indep-pairwise 2000 50 0.1 --out indep

wc -l indep.prune.in
plink2 --bfile data_maf0.01_rs_ref --extract indep.prune.in --out data_maf0.01_rs_ref_pruned --make-just-bim


## Prune hm3 from chromosome X

# Get 1kg vcf
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# Create psam
sed 1d integrated_call_samples_v3.20130502.ALL.panel | awk '{ print $1, $1, "0", "0", $4, "0" }' | sed 's/female/2/g' | sed 's/male/1/g'  > 1kg_chrX.psam
head 1kg_chrX.psam

# check same number of samples
wc -l 1kg_chrX.psam
zcat ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz | grep -m 1 "#CHROM" | tr "\t" "\n" | grep -E [0-9] > vcf_idlist
wc -l vcf_idlist

# Get hm3 positions etc from chromosome X
grep -wf hapmap3_r1_b36_fwd_consensus.qc.poly.recode.snp data.chrX.snp-stats | awk '{ print $3, $4, $4, $2 }' > hm3_x.bed

# Extract from vcf and prune
plink2 --vcf ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz --psam 1kg_chrX.psam --extract range hm3_x.bed --set-all-var-ids @:#:[b37]\$r,\$a --new-id-max-allele-len 161 --max-alleles 2 --indep-pairwise 2000 50 0.1 --out indep_x --split-par hg19

# Make bed
cut -d ":" -f 2 indep_x.prune.in > hm3_x_prune_pos
grep -wf hm3_x_prune_pos hm3_x.bed | awk '{ print "chr"$1, $2, $3, $4 }' > hm3_x_prune.bed


# Combine with autosomal prune list
awk '{ print "chr"$1, $4, $4, $2 }' data_maf0.01_rs_ref_pruned.bim > hm3_auto_prune_b37.bed

cat hm3_auto_prune_b37.bed hm3_x_prune.bed > hm3_prune_b37.bed

# Add tophits
> tophitsnps.txt
for f in ../../resources/genotypes/tophits/*.txt
do
    awk '{ print $1 }' $f | cut -d "_" -f 1 ${prunefile} | tr ":" " " | awk '{ print "chr"$1, $2, $2, $1":"$2 }' >> tophitsnps.txt
done
cat tophitsnps.txt hm3_prune_b37.bed | sort | uniq > hm3_prune_th_b37.bed
wc -l hm3_prune_th_b37.bed
cut -d " " -f 1 hm3_prune_th_b37.bed | sort | uniq -c

# Liftover to hg38
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
chmod 755 liftOver
gunzip hg19ToHg38.over.chain.gz

./liftOver hm3_prune_th_b37.bed hg19ToHg38.over.chain hm3_prune_th_b38.bed unmapped

wc -l unmapped

sed 's/chr//g' hm3_prune_th_b37.bed > hm3_prune_th_b37_nochr.bed
sed 's/chr//g' hm3_prune_th_b38.bed > hm3_prune_th_b38_nochr.bed

gzip -c hm3_prune_th_b37_nochr.bed > ../../resources/genotypes/hm3_prune_th_hg19.bed.gz
gzip -c hm3_prune_th_b38_nochr.bed > ../../resources/genotypes/hm3_prune_th_hg38.bed.gz

head hm3_prune_th_b38_nochr.bed
head hm3_prune_th_b37_nochr.bed



# wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz
# wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz.md5
# wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz.tbi


