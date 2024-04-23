library(data.table)
# library(fst)

fn <- commandArgs(T)[1]
stopifnot(file.exists(fn))
a <- data.table::fread(fn, header=TRUE)

# Harmonise alleles to alphabetical

ord <- a$A1 > a$A2
table(ord)
a$BETA[ord] <- a$BETA[ord] * -1
a$AF1[ord] <- 1-a$AF1[ord]

b <- a[, .(SNP, BETA, SE, AF1, N)]
head(b)
# p <- paste0(fn, ".fst")
# fst::write_fst(a, path=p, compress=100)
p <- paste0(fn, ".rds")
saveRDS(b, p)

## Note
# Was previously using fst but finding it can be cumbersome to install
# So reverting to using rds as it

## redundant - these are bigger than compressed files
# writebingwas <- function(a, fn) {
#     con <- file(fn, "wb")
#     n <- nrow(a)
#     writeBin(n, con)
#     writeBin(a[, BETA], con)
#     writeBin(a[, SE], con)
#     writeBin(a[, AF1], con)
#     writeBin(a[, N], con)
#     close(con)
# }

# readbingwas <- function(fn) {
#     con <- file(fn, "rb")
#     n <- readBin(con, integer(), n = 1)
#     tibble(
#         BETA = readBin(con, numeric(), n=n),
#         SE = readBin(con, numeric(), n=n),
#         AF1 = readBin(con, numeric(), n=n),
#         N = readBin(con, integer(), n=n)
#     )
#     close(con)
# }
# writebingwas(a, "temp.bin")
