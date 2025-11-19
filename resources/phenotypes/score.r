library(dplyr)
library(data.table)
nthreads <- as.numeric(Sys.getenv("env_threads"))


args <- commandArgs(T)
phenolist <- scan(args[1], what="character")
score_dir <- args[2]
out_dir <- args[3]

phenolist1 <- phenolist
phenolist <- gsub("bioavail_testosterone", "bioavail-testosterone", phenolist)
phenolist <- gsub("bioavail_", "bioavail-testosterone_", phenolist)
dat <- tibble(
    fn=phenolist,
    bn=basename(fn) %>% gsub(".phen$", "", .),
    phen=sapply(bn, \(x) strsplit(x, "_")[[1]][1]),
    agebin=sapply(bn, \(x) strsplit(x, "_")[[1]][2]),
    r = NA,
    n = NA,
    b = NA,
    se = NA,
    pval = NA,
    phen_sd = NA,
    phen_m = NA,
    score_sd = NA,
    score_m = NA,
)

dat$phen <- gsub("bioavail-testosterone", "bioavail_testosterone", dat$phen)
dat$fn <- phenolist1

print(unique(dat$phen))

for(i in 1:nrow(dat)) {
    a <- fread(dat$fn[i], nThread = nthreads, keepLeadingZeros=TRUE)
    a <- a[,1:3]
    names(a)[1:3] <- c("FID", "IID", "value")
    cf <- gsub(".phen$", ".covs", dat$fn[i])
    covs <- fread(cf, nThread = nthreads, keepLeadingZeros=TRUE)
    a <- inner_join(a, covs, by=c("FID"="V1", "IID"="V2"))
    b <- fread(file.path(score_dir, paste0(dat$phen[i], ".sscore")), nThread = nthreads, keepLeadingZeros=TRUE)
    a$FID <- as.character(a$FID)
    b$`#FID` <- as.character(b$`#FID`)
    ab <- inner_join(a, b, by=c("FID"="#FID", "IID"="IID"))
    f <- paste0("value ~ SCORE1_AVG + ", paste(names(covs)[-c(1:2)], collapse = " + "))
    r <- cor(ab$value, ab$SCORE1_AVG)
    m <- summary(lm(f, ab))
    dat$n[i] <- m$df[1] + m$df[2]
    dat$r[i] <- r
    dat$b[i] <- m$coef[2,1]
    dat$se[i] <- m$coef[2,2]
    dat$pval[i] <- m$coef[2,4]
    dat$phen_sd[i] <- sd(ab$value, na.rm=T)
    dat$phen_m[i] <- mean(ab$value, na.rm=T)
    dat$score_sd[i] <- sd(ab$SCORE1_AVG, na.rm=T)
    dat$score_m[i] <- mean(ab$SCORE1_AVG, na.rm=T)
}

dat <- dat %>% dplyr::select(!fn)
print(head(dat))

saveRDS(dat, file=file.path(out_dir, "scores.rds"))
