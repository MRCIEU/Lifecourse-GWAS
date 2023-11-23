library(data.table)
library(dplyr)
library(here)
library(glue)

# read in bim file

args <- commandArgs(T)
bfile <- args[1]
outfile <- args[2]

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
    "{here('bin', 'plink2')} --bfile {bfile} --ref-allele {switchfile} 2 1 --make-bed --out {outfile} --keep-allele-order"
) %>% system()

# make copy of original bim file

file.copy(from=paste0(outfile, ".bim"), to=paste0(outfile, ".bim.orig"))

# update variant IDs

bim <- data.table::fread(paste0(outfile, ".bim"))
ch <- bim$V5 < bim$V6
table(ch)
stopifnot(all(ch))

bim$V2 <- paste0(bim$V1, ":", bim$V4, "_", bim$V5, "_", bim$V6)

data.table::fwrite(bim, file=paste0(outfile, ".bim"), quote=FALSE, col.names=FALSE, sep=" ")
