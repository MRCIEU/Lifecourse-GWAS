library(data.table)
library(dplyr)
library(dotenv)
library(digest)
library(glue)

get_lambda <- function(pvector, pl) {
    pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
    o <- -log10(sort(pvector, decreasing=FALSE))
    e <- -log10(ppoints(length(pvector)))
    cs <- qchisq(1-pvector, 1)
    lambda <- median(cs, na.rm=TRUE) / qchisq(0.5, 1)
    png(pl)
    plot(e, o, xlab="Expected -log10(p)", ylab="Observed -log10(p)", pch=20, col="blue", main=paste("null lambda =", lambda, "; nsnp = ", length(pvector)))
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

selvariants <- fread(aname, header=TRUE) %>% as_tibble() 

str(selvariants)

min(selvariants$P, na.rm=T)

get_lambda(selvariants$P, file.path(resdir, "null_qq.png"))
