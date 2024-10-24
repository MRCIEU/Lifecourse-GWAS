library(dplyr)
library(data.table)
library(parallel)


get_lambda <- function(pvector, pl=FALSE) {
    pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
    o <- -log10(sort(pvector, decreasing=FALSE))
    e <- -log10(ppoints(length(pvector)))
    if(pl) {
        plot(e, o, xlab="Expected -log10(p)", ylab="Observed -log10(p)", pch=20, col="blue")
        abline(0, 1, col="red")
    }
    cs <- qchisq(1-pvector, 1)
    lambda <- median(cs, na.rm=TRUE) / qchisq(0.5, 1)
    return(lambda)
}

parse_file <- function(gwas, resdir) {
    message("Reading data")
    a <- lapply(gwas$fn, data.table::fread, header=TRUE, data.table=FALSE, nThread=1) %>% bind_rows() %>% as_tibble()
    # a <- fread(fn, header=TRUE, data.table=FALSE, colClasses=c()) %>% as_tibble()
    # a[,2] <- as.integer(a[,2, drop=TRUE])
    # a[,6] <- as.numeric(a[,6, drop=TRUE])
    # a[,7] <- as.numeric(a[,7, drop=TRUE])
    # a[,8] <- as.numeric(a[,8, drop=TRUE])
    # a[,10] <- as.numeric(a[,10, drop=TRUE])
    # a[,11] <- as.numeric(a[,11, drop=TRUE])
    # a[,12] <- as.numeric(a[,12, drop=TRUE])
    # a[,13] <- as.numeric(a[,13, drop=TRUE])

    message("Aligning")
    a <- subset(a, !is.na(a$GENPOS))
    flip <- a$ALLELE1 > a$ALLELE0
    temp <- a$ALLELE1[flip]
    a$ALLELE1[flip] <- a$ALLELE0[flip]
    a$ALLELE0[flip] <- temp
    a$BETA[flip] <- a$BETA[flip] * -1
    a$A1FREQ[flip] <- 1 - a$A1FREQ[flip]
    
    message("IDs")
    a$ID <- paste0(a$CHROM, ":", a$GENPOS, "_", a$ALLELE1, "_", a$ALLELE0)
    message("P-values")
    a$pval <- 10^-a$LOG10P

    message("Summarising")
    nom <- gwas$gwas[1] %>% strsplit("_") %>% unlist
    s <- tibble(
        phen = nom[1],
        age = nom[2],
        sex = nom[3],
        nsnp = nrow(a),
        meanchisq = mean(a$CHISQ, na.rm=TRUE),
        medianchisq = median(a$CHISQ, na.rm=TRUE),
        lambda1 = medianchisq / qchisq(0.5, 1),
        # lambda2 = get_lambda(a$pval),
        minp = min(a$pval, na.rm=TRUE),
        max_n = max(a$N, na.rm=TRUE),
        min_n = min(a$N, na.rm=TRUE),
        mean_n = mean(a$N, na.rm=TRUE)
    )
    message("Writing")
    b <- a %>% dplyr::select(chr=CHROM, pos=GENPOS, id=ID, ea=ALLELE1, oa=ALLELE0, eaf=A1FREQ, beta=BETA, se=SE, n=N)
    out <- file.path(resdir, paste0(gwas$gwas[1], ".regenie.gz"))
    data.table::fwrite(b, file=out, row=FALSE, col=TRUE, qu=FALSE, sep="\t", nThread=1)

    return(s)
}

args <- commandArgs(T)

resdir <- args[1]
nchr <- as.numeric(args[2])
phenolist <- args[3]

# resdir <- "/local-scratch/projects/Lifecourse-GWAS/gib/alspac/results3/04"
# nchr <- 22
# phenolist <- scan("/local-scratch/projects/Lifecourse-GWAS/gib/alspac/phen_proc2/phenolist", what="character")
phenolist <- scan(phenolist, what="character")

gwass <- lapply(phenolist, \(x) {
    nom <- basename(x) %>% gsub(".phen$", "", .)
    d <- dirname(x)
    tibble(
        gwas = nom,
        fn = file.path(d, "regenie", paste0("step2_", 1:nchr, "_", nom, ".regenie.gz")),
        exists = file.exists(fn)
    )
}) %>% bind_rows()
table(gwass$exists)

if(!all(gwass$exists)) {
    stop("The following files do not exist: ", paste0(gwass$fn[!gwass$exists], collapse="\n"))
}


ugwas <- unique(gwass$gwas)
gwas_summary <- mclapply(ugwas, \(x) {
    message(x)
    parse_file(subset(gwass, gwas == x), resdir)
}, mc.cores=as.numeric(Sys.getenv("env_threads")))

saveRDS(gwas_summary, file=file.path(resdir, "gwas_summary.rds"))
