library(data.table)
library(dplyr)
library(dotenv)

standardise <- function(d, ea_col="ea", oa_col="oa", beta_col="beta", eaf_col="eaf", chr_col="chr", pos_col="pos", vid_col="vid")
{
    toflip <- d[[ea_col]] > d[[oa_col]]
    # d[[eaf_col]][toflip] <- 1 - d[[eaf_col]][toflip]
    # d[[beta_col]][toflip] <- d[[beta_col]][toflip] * -1
    temp <- d[[oa_col]][toflip]
    d[[oa_col]][toflip] <- d[[ea_col]][toflip]
    d[[ea_col]][toflip] <- temp
    d[[vid_col]] <- paste0(d[[chr_col]], ":", d[[pos_col]], "_", d[[ea_col]], "_", d[[oa_col]])
    d
}

load_dot_env("config.env")
Sys.getenv()

cohort <- Sys.getenv("cohort_name")
minmaf <- as.numeric(Sys.getenv("env_minmaf"))
mininfo <- as.numeric(Sys.getenv("env_mininfo"))

args <- commandArgs(T)

aname <- "fastgwa_bgen_2.fastGWA"
aname <- args[1]

resdir <- "/local-scratch/projects/Lifecourse-GWAS/gib/alspac/results6/00"
resdir <- args[2]

inclist <- "/local-scratch/projects/Lifecourse-GWAS/gib/alspac/geno_proc6/variant_inclusion.txt"
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
table(selvariants$CHR)

min(selvariants$P, na.rm=T)

selvariants <- standardise(selvariants, ea_col="A1", oa_col="A2", beta_col="BETA", eaf_col="AF1", chr_col="CHR", pos_col="POS", vid_col="VID")
table(duplicated(selvariants$SNP))
table(duplicated(selvariants$VID))

head(selvariants)

keepvars <- unique(selvariants$SNP)

selvariants <- dplyr::select(selvariants, CHR, POS, VID, SNP, A1, A2, AF1, INFO)
fwrite(selvariants, file=file.path(resdir, "variants.txt.gz"), row.names=F, col.names=T, quote=F, sep="\t")

write.table(keepvars, file=inclist, row.names=F, col.names=F, quote=F)
