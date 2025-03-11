library(data.table)
library(dplyr)
library(R.utils)
library(dotenv)
load_dot_env("config.env")

create_variantid <-function(chr,pos,a1,a2) {
  toflip <- a1 > a2
  temp <- a1[toflip]
  a1[toflip] <- a2[toflip]
  a2[toflip] <- temp

  # create hashes when alleles nchar > 10
  # allele ea
  nch1 <- nchar(a1) > 10
  nch2 <- nchar(a2) > 10
  if(any(nch1)) {
    index <- which(nch1)
    a1[nch1] <- sapply(a1[nch1], \(allele) {
        digest::digest(allele, algo="murmur32")
    })
  }

  if(any(nch2)) {
    index <- which(nch2)
    a2[nch2] <- sapply(a2[nch2], \(allele) {
        digest::digest(allele, algo="murmur32")
    })
  }

  #  create variantid
  variantid <- paste0(chr, ":", pos, "_", a1, "_", a2)

  return(variantid)
}

# create_variantid(
#     c(1,2,3),
#     c(10000, 10000, 10000),
#     c("G", "TACTGTGTGTGTGTGTGT", "AAAAA"),
#     c("A", "GACTGTGTGTGTGTGTGT", "AAAAAG")
# )

standardise <- function(d, ea_col="ea", oa_col="oa", beta_col="beta", eaf_col="eaf", chr_col="chr", pos_col="pos", vid_col="vid")
{
    toflip <- d[[ea_col]] > d[[oa_col]]
    d[[eaf_col]][toflip] <- 1 - d[[eaf_col]][toflip]
    d[[beta_col]][toflip] <- d[[beta_col]][toflip] * -1
    a1 <- d[[oa_col]]
    a2 <- d[[ea_col]]
    temp <- d[[oa_col]][toflip]
    d[[oa_col]][toflip] <- d[[ea_col]][toflip]
    d[[ea_col]][toflip] <- temp

    d[[vid_col]] <- create_variantid(d[[chr_col]], d[[pos_col]], d[[ea_col]], d[[oa_col]])
    d
}

get_lambda <- function(pvector) {
    pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
    o <- -log10(sort(pvector, decreasing=FALSE))
    e <- -log10(ppoints(length(pvector)))
    cs <- qchisq(1-pvector, 1)
    lambda <- median(cs, na.rm=TRUE) / qchisq(0.5, 1)
    return(lambda)
}

process_gwas <- function(fn, a) {
    fn2 <- basename(fn) %>% gsub(".fastGWA", "", .)
    nom <- fn2 %>% strsplit("_") %>% unlist
    
    s <- tibble(
        phen = nom[1],
        age = nom[2],
        sex = nom[3],
        nsnp = nrow(a),
        meanpval = mean(a$P, na.rm=TRUE),
        medianpval = median(a$P, na.rm=TRUE),
        lambda = qchisq(medianpval, df=1, low=FALSE) / qchisq(0.5, 1),
        minp = min(a$P, na.rm=TRUE),
        max_n = max(a$N, na.rm=TRUE),
        min_n = min(a$N, na.rm=TRUE),
        mean_n = mean(a$N, na.rm=TRUE)
    )
    return(s)
}

fn <- commandArgs(T)[1]
ref <- commandArgs(T)[2]
nthreads <- as.numeric(Sys.getenv("env_threads"))
stopifnot(file.exists(fn))

message("reading gwas")
a <- data.table::fread(fn, header=TRUE, nThread = nthreads) %>% as_tibble()
a <- subset(a, !is.na(P) & !is.nan(P) & is.finite(P) & P <= 1 & P >= 0)

if(Sys.getenv("genome_build") == "hg19") {
    stopifnot(file.exists(ref))
    message("reading reference")
    v <- data.table::fread(ref, header=TRUE, nThread = nthreads) %>% as_tibble()
    names(a)[names(a) == "POS"] <- "POS19"
    a <- inner_join(v, a)
}

message("harmonising")
a <- standardise(a, ea_col="A1", oa_col="A2", beta_col="BETA", eaf_col="AF1", chr_col="CHR", pos_col="POS", vid_col="VID")
s <- process_gwas(fn, a)
a <- a %>% dplyr::select(VID, BETA, SE, EAF=AF1, N, P)
fwrite(a, file=paste0(fn, ".gz"))
saveRDS(s, file=paste0(fn, ".summary.rds"))
