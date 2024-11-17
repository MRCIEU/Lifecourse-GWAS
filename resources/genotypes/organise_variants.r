library(data.table)
library(dplyr)
library(dotenv)
library(digest)
library(glue)

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
    # d[[beta_col]][toflip] <- d[[beta_col]][toflip] * -1
    a1 <- d[[oa_col]]
    a2 <- d[[ea_col]]
    temp <- d[[oa_col]][toflip]
    d[[oa_col]][toflip] <- d[[ea_col]][toflip]
    d[[ea_col]][toflip] <- temp

    d[[vid_col]] <- create_variantid(d[[chr_col]], d[[pos_col]], d[[ea_col]], d[[oa_col]])
    d
}

get_lambda <- function(pvector, pl) {
    pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
    o <- -log10(sort(pvector, decreasing=FALSE))
    e <- -log10(ppoints(length(pvector)))
    cs <- qchisq(1-pvector, 1)
    lambda <- median(cs, na.rm=TRUE) / qchisq(0.5, 1)
    png(pl)
    plot(e, o, xlab="Expected -log10(p)", ylab="Observed -log10(p)", pch=20, col="blue", main=paste("random phenotype lambda =", lambda))
    abline(0, 1, col="red")
    dev.off()
    return(lambda)
}


load_dot_env("config.env")

cohort <- Sys.getenv("cohort_name")
minmaf <- as.numeric(Sys.getenv("env_minmaf"))
mininfo <- as.numeric(Sys.getenv("env_mininfo"))

args <- commandArgs(T)
aname <- args[1]
resdir <- args[2]
inclist <- args[3]

allvariants <- fread(aname, header=TRUE) %>% as_tibble() 

# plot info distribution
# plot af distribution
# filter by maf / info
# generate new variant id
# write keep file
# write variant index file - chr, pos, newid, oldif, a1, a2, af, info



png(file.path(resdir, "info_dist.png"))
hist(allvariants$INFO, breaks=100, main=paste(cohort, "INFO", nrow(allvariants), "variants"), xlab="INFO")
dev.off()

png(file.path(resdir, "af_dist.png"))
hist(allvariants$AF1, breaks=100, main=paste(cohort, "AF", nrow(allvariants), "variants"), xlab="INFO")
dev.off()

selvariants <- subset(allvariants, INFO > mininfo & AF1 > minmaf & AF1 < 1-minmaf)
dim(selvariants)
table(allvariants$CHR) %>% as.data.frame
table(selvariants$CHR) %>% as.data.frame

min(selvariants$P, na.rm=T)

get_lambda(selvariants$P, file.path(resdir, "null_qq.png"))

if(Sys.getenv("genome_build") == "hg19") {
  names(selvariants)[names(selvariants) == "POS"] <- "POS19"
  system("gunzip -c bin/liftOver.gz > bin/liftOver")
  system("chmod 755 bin/liftOver")
  a <- dplyr::select(selvariants, CHR, POS1=POS19, POS2=POS19, SNP) %>% mutate(CHR = paste0("chr", CHR))
  tf <- tempfile()
  fwrite(a, file=tf, row=F, col=F, qu=F, sep="\t")
  cmd <- glue("./bin/liftOver {tf} resources/genotypes/hg19ToHg38.over.chain {tf}.out {tf}.unlifted")
  system(cmd)
  b <- fread(glue("{tf}.out"), he=F)
  b <- b %>% filter(V1 %in% paste0("chr", c(1:22, "X"))) %>% mutate(V1 = gsub("chr", "", V1))
  table(b$V1) %>% as.data.frame()
  ab <- left_join(b, a, by=c("V4"="SNP"))
  ab$CHR <- gsub("chr", "", ab$CHR)
  ab$CHR[ab$CHR == 23] <- "X"

  # Remove variants that don't match the same chromosome
  ab <- subset(ab, ab$CHR == ab$V1)
  ab <- ab %>% dplyr::select(SNP=V4, POS=V2)
  ab <- subset(ab, !duplicated(SNP))
  selvariants <- left_join(ab, selvariants, by="SNP")
  selvariants
}

selvariants <- standardise(selvariants, ea_col="A1", oa_col="A2", beta_col="BETA", eaf_col="AF1", chr_col="CHR", pos_col="POS", vid_col="VID")
selvariants <- selvariants %>%
  arrange(desc(INFO)) %>%
  filter(!duplicated(VID)) %>%
  arrange(CHR, POS)

table(duplicated(selvariants$SNP))
table(duplicated(selvariants$VID))
head(selvariants)

keepvars <- unique(selvariants$SNP)

if(Sys.getenv("genome_build") == "hg19") {
  selvariants <- dplyr::select(selvariants, CHR, POS, POS19, VID, SNP, A1, A2, AF1, INFO)
} else {
  selvariants <- dplyr::select(selvariants, CHR, POS, VID, SNP, A1, A2, AF1, INFO)
}
fwrite(selvariants, file=file.path(resdir, "variants.txt.gz"), row.names=F, col.names=T, quote=F, sep="\t")

write.table(keepvars, file=inclist, row.names=F, col.names=F, quote=F)
