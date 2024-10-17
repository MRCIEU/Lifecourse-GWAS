library(data.table)
library(dplyr)
library(R.utils)

# fn <- "/local-scratch/projects/Lifecourse-GWAS/gib/alspac/results/03/bmi_10-11.fastGWA"
fn <- commandArgs(T)[1]
ref <- commandArgs(T)[2]
rem <- commandArgs(T)[3]
nthread <- as.numeric(commandArgs(T)[4])
stopifnot(file.exists(fn))
stopifnot(file.exists(ref))
stopifnot(file.exists(rem))

message("reading gwas")
a <- data.table::fread(fn, header=TRUE, nThreads = nthread) %>% as_tibble()
message("reading reference")
v <- data.table::fread(ref, header=TRUE, nThreads = nthread) %>% as_tibble()

# Sometimes GCTA doesn't remove all the variants that it's supposed to
# Do this manually
rem <- scan(rem, "character")
remflag <- sum(rem %in% a$SNP)
if(remflag > 0) {
    message(sum(remflag), " flagged variants need to be removed from GWAS")
}
a <- subset(a, !SNP %in% rem)

check_against_reference <- function(v, a) {
    if(all(a$SNP == v$SNP)) {
        message("gwas matches reference")
        return(a)
    } else {
        message("aligning gwas with reference")
        print(dim(a))
        a <- subset(a, SNP %in% v$SNP)
        print(dim(a))
        dum <- subset(v, !SNP %in% a$SNP)
        print(dum)
        dum$BETA <- NA_real_
        dum$SE <- NA_real_
        dum$P <- NA_real_
        dum$N <- NA_real_
        dum$CHR <- as.character(dum$CHR)
        a$CHR <- as.character(a$CHR)
        dum <- tibble(CHR=as.character(dum$CHR), SNP=dum$SNP, A2=dum$EA, A1=dum$OA, AF1=dum$EAF)
        print(dum)
        a <- bind_rows(a, dum)
        m <- match(v$SNP, a$SNP)
        a <- a[m,]
        stopifnot(all(a$SNP == v$SNP))
        return(a)
    }
}

a <- check_against_reference(v, a)

# Harmonise alleles to alphabetical
message("harmonising")
ord <- a$A1 > a$A2
table(ord)
a$BETA[ord] <- a$BETA[ord] * -1
a$AF1[ord] <- 1-a$AF1[ord]
temp <- a$A1[ord]
a$A1[ord] <- a$A2[ord]
a$A2[ord] <- temp

message("writing")
b <- a %>% dplyr::select(BETA, SE)
saveRDS(b, file=paste0(fn, ".rds"))
