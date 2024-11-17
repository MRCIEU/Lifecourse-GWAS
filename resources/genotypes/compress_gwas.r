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


# fn <- "/local-scratch/projects/Lifecourse-GWAS/gib/alspac/results/03/bmi_10-11.fastGWA"
fn <- commandArgs(T)[1]
ref <- commandArgs(T)[2]
nthreads <- as.numeric(Sys.getenv("env_threads"))
stopifnot(file.exists(fn))
stopifnot(file.exists(ref))

message("reading gwas")
a <- data.table::fread(fn, header=TRUE, nThread = nthreads) %>% as_tibble()

if(Sys.getenv("genome_build") == "hg19") {
    message("reading reference")
    v <- data.table::fread(ref, header=TRUE, nThread = nthreads) %>% as_tibble()
    names(a)[names(a) == "POS"] <- "POS19"
    a <- left_join(v, a)
}

message("harmonising")
a <- standardise(a, ea_col="A1", oa_col="A2", beta_col="BETA", eaf_col="AF1", chr_col="CHR", pos_col="POS", vid_col="VID")
a <- a %>% dplyr::select(VID, BETA, SE, EAF=AF1, N, P)
fwrite(a, file=paste0(fn, ".gz"))
