library(microbenchmark)
library(data.table)
library(dplyr)
library(fst)

fn <- "/local-scratch/projects/Lifecourse-GWAS/gib/alspac/results/03/bmi_10-11.fastGWA"
a <- data.table::fread(fn, header=TRUE) %>% as_tibble()
v <- a %>% select(CHR, )
ss <- b %>% 
ssmin <- bs
ssmin


x <- expand.grid(
    what = c("variants", "ss", "ssmin"),
    format = c("rds", "data.table gzip", "data.table txt", "fst"),
    file_size = NA,
    read_time = NA,
    write_time = NA
)

i <- 1
for(i in which(x$what == "ss")) {
    message(i)
    fn <- tempfile()
    if(x$format[i] == "rds") {
        message("rds")
        obj <- get(as.character(x$what[i]))
        x$write_time[i] <- microbenchmark(saveRDS(obj, file=fn), times=5) %>% summary %>% as.data.frame %>% {.$median}
        x$file_size[i] <- file.size(fn)
        x$read_time[i] <- microbenchmark(temp <- readRDS(file=fn), times=5) %>% summary %>% as.data.frame %>% {.$median}
    }
    if(x$format[i] == "data.table txt") {
        message("data.table txt")
        obj <- get(as.character(x$what[i]))
        x$write_time[i] <- microbenchmark(fwrite(obj, file=fn, compress="none"), times=5) %>% summary %>% as.data.frame %>% {.$median}
        x$file_size[i] <- file.size(fn)
        x$read_time[i] <- microbenchmark(temp <- fread(file=fn), times=5) %>% summary %>% as.data.frame %>% {.$median}
    }
    if(x$format[i] == "data.table gzip") {
        message("data.table gzip")
        fn2 <- paste0(fn, ".gz")
        obj <- get(as.character(x$what[i]))
        x$write_time[i] <- microbenchmark(fwrite(obj, file=fn2, compress="gzip"), times=5) %>% summary %>% as.data.frame %>% {.$median}
        x$file_size[i] <- file.size(fn2)
        x$read_time[i] <- microbenchmark(temp <- fread(file=fn2), times=5) %>% summary %>% as.data.frame %>% {.$median}
    }
    if(x$format[i] == "fst") {
        message("data.table fst")
        obj <- get(as.character(x$what[i]))
        x$write_time[i] <- microbenchmark(write_fst(obj, path=fn, compress=100), times=5) %>% summary %>% as.data.frame %>% {.$median}
        x$file_size[i] <- file.size(fn)
        x$read_time[i] <- microbenchmark(temp <- read_fst(path=fn), times=5) %>% summary %>% as.data.frame %>% {.$median}
    }

}





# Summary
# - fst vs rds
#     is twice as fast at writing and reading, but only tested on SSD. 3 seconds vs 6 seconds
#     - it might be tricky to install for some
#     - file sizes are about 10-20% larger
# - fread txt for variants is comparable to rds, worth using for variants so that it's easy to read




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
