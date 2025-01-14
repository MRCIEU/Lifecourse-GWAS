library(data.table)
library(dplyr)
library(here)
library(digest)
nthreads <- as.numeric(Sys.getenv("env_threads"))

# read in bim file
args <- commandArgs(T)
bfile <- args[1]
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
file.copy(paste0(bfile, ".bim"), paste0(bfile, ".bim.orig"))
bim <- data.table::fread(paste0(bfile, ".bim.orig"), nThread = nthreads)
bim$switch <- bim$V5 > bim$V6
temp <- bim$V5[bim$switch]
bim$A1 <- bim$V5
bim$A2 <- bim$V6
bim$A1[bim$switch] <- bim$A2[bim$switch]
bim$A2[bim$switch] <- temp
table(bim$switch)
str(bim)

# update variant IDs
ch <- bim$A1 < bim$A2
table(ch)
stopifnot(all(ch))

bim$V2 <- paste0(bim$V1, ":", bim$V4, "_", compress_alleles(bim$A1), "_", compress_alleles(bim$A2))
bim <- subset(bim, select=-c(A1, A2, switch))

# manage duplicates

# from https://github.com/bhklab/genefu
rename.duplicate <- function (x, sep="_", verbose=FALSE) {
	x <- as.character(x)
	duplix <- duplicated(x)
	duplin <- x[duplix]

	ix <- numeric(length=length(unique(duplin)))
	names(ix) <- unique(duplin)
	retval <- numeric(length=length(duplin))
	for(i in 1:length(duplin)) { retval[i] <- ix[duplin[i]] <- ix[duplin[i]] + 1 }
	retval <- retval + 1
	x[duplix] <- paste(duplin, retval, sep=sep)

	if (verbose) { message(sprintf("%i duplicated names", length(duplin))) }
	
	return (x)
}

bim$V2 <- rename.duplicate(bim$V2, "_duplicate")

data.table::fwrite(bim, file=paste0(bfile, ".bim"), quote=FALSE, col.names=FALSE, sep=" ", nThread = nthreads)


