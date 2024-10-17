library(dplyr)
library(data.table)

args <- commandArgs(T)
rem <- scan(args[1], what="character")
out <- args[2]
nthreads <- as.numeric(Sys.getenv("env_threads"))

afreqlist <- scan(args[3], what="character")

l <- lapply(afreqlist, \(fn) {
    a <- fread(fn, nThread = nthreads) %>% as_tibble()
    names(a) <- c("CHR", "SNP", "OA", "EA", "PROVREF", "EAF", "OBS_CT")
    a <- subset(a, !SNP %in% rem)
    stopifnot(all(a$PROVREF == "Y"))
    a <- subset(a, select=-c(PROVREF))
    a$flipped <- a$EA > a$OA
    table(a$flipped)
    a$EAF[a$flipped] <- 1 - a$EAF[a$flipped]
    temp <- a$EA[a$flipped]
    a$EA[a$flipped] <- a$OA[a$flipped]
    a$OA[a$flipped] <- temp
    a$CHR <- as.character(a$CHR)
    str(a)
    a
}) %>% bind_rows()

fwrite(l, out)
