library(data.table)
library(dplyr)
library(here)
library(glue)
library(digest)

# read in bim file

args <- commandArgs(T)
bfile <- args[1]
outfile <- args[2]
plinkbin <- args[3]

compress_alleles <- function(a) {
    i <- nchar(a) > 10
    if(any(i)) {
        a[i] <- sapply(a[i], \(x) digest(x, algo="murmur32"))
    }
    return(a)
}

# bim <- tribble(
#     ~V1, ~V2, ~V3, ~V4, ~V5, ~V6,
#     1, "rs234", 0, 1000, "G", "A",
#     1, "rs234", 0, 1000, "GGGGGGGGGGGCCA", "A",
#     1, "rs235", 0, 2000, "GGGGGGGGGGGCCAG", "A",
#     1, "rs237", 0, 3000, "G", "T"
# )

# update orig bim
bim <- data.table::fread(paste0(bfile, ".bim"))
bim$switch <- bim$V5 > bim$V6
temp <- bim$V5[bim$switch]
bim$V5[bim$switch] <- bim$V6[bim$switch]
bim$V6[bim$switch] <- temp
table(bim$switch)

# switch alleles in data

switchfile <- paste0(outfile, ".switch")
data.table::fwrite(subset(bim, select=c(V2, V6)), file=switchfile, quote=FALSE, col.names=FALSE, sep=" ")

glue(
    "{plinkbin} --bfile {bfile} --rm-dup force-first --ref-allele {switchfile} 2 1 --make-bed --out {outfile} --keep-allele-order"
) %>% system()

# make copy of original bim file

file.copy(from=paste0(outfile, ".bim"), to=paste0(outfile, ".bim.orig"))

# update variant IDs

bim <- data.table::fread(paste0(outfile, ".bim"))
ch <- bim$V5 < bim$V6
table(ch)
stopifnot(all(ch))

bim$V2 <- paste0(bim$V1, ":", bim$V4, "_", compress_alleles(bim$V5), "_", compress_alleles(bim$V6))

data.table::fwrite(bim, file=paste0(outfile, ".bim"), quote=FALSE, col.names=FALSE, sep=" ")
