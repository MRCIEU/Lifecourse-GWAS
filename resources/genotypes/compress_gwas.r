library(data.table)
library(fst)

fn <- commandArgs(T)[1]
stopifnot(file.exists(fn))
a <- data.table::fread(fn, header=TRUE)
names(a) <- c("SNP", "BETA", "SE", "AF1", "N")
b <- a[, .(SNP, BETA, SE, AF1, N)]
b
p <- paste0(fn, ".fst")
fst::write_fst(a, path=p, compress=100)


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
