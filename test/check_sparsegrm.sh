## code to keep only a subset for testing purposes
# head -n 1000 data_chr01.fam > temp1
# tail -n 1000 data_chr01.fam > temp2
# cat temp1 temp2 | awk '{print $1, $2}' > keep

# for f in *bed
# do
#     echo $f
#     bn=$(echo $f | sed 's/.bed//g')
#     plink2 --bfile $bn --keep keep --make-bed --out subset/${bn}
# done


# This code generates the sparse GRM the long way - first full grm then making it sparse
# expect it gives similar answer as the quick version in the code
bin/gcta64 --bfile ${genotype_processed_dir}/scratch/indep --make-grm-bin --out ${genotype_processed_dir}/${bfile_prefix}_test
bin/gcta64 --grm ${genotype_processed_dir}/${bfile_prefix}_test --make-bK-sparse 0.05 --out ${genotype_processed_dir}/${bfile_prefix}_test

# check in R:
# a <- read.table("/mnt/storage/private/mrcieu/research/scratch/Lifecourse-GWAS/gib/geno_proc/data_chr_test.grm.sp")
# b <- read.table("/mnt/storage/private/mrcieu/research/scratch/Lifecourse-GWAS/gib/geno_proc/data_chr.grm.sp")

# cor(a$V3, b$V3)

# ab <- merge(a, b, by.x=c("V1", "V2"), by.y=c("V1", "V2"))

# dim(ab)

# head(ab)
# cor(ab$V3.x, ab$V3.y)
# ab2 <- subset(ab, V3.y != 1)
# cor(ab2$V3.x, ab2$V3.y)
# ab2
