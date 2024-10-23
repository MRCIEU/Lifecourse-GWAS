library(dplyr)
library(data.table)

a <- fread("regenie_test_bgen_bmi_40-45_both.regenie", header = TRUE, data.table = FALSE) %>% as_tibble()

b <- fread("regenie_test_bestguess_bmi_40-45_both.regenie", header = TRUE, data.table = FALSE) %>% as_tibble()

table(duplicated(paste(a$CHROM, a$GENPOS)))
table(duplicated(paste(b$CHROM, b$GENPOS)))




a$chrpos <- paste(a$CHROM, a$GENPOS)
b$chrpos <- paste(b$CHROM, b$GENPOS)

ab <- inner_join(a, b, by="chrpos")

ab <- subset(ab, ALLELE0.x == ALLELE1.y & ALLELE1.x == ALLELE0.y)


cor(ab$BETA.x, ab$BETA.y)

cor(ab$SE.x, ab$SE.y)

summary(lm(BETA.x ~ BETA.y, data=ab))

summary(lm(SE.x ~ SE.y, data=ab))




png("../Lifecourse-GWAS-alspac/beta_bg_bgen.png")
plot(BETA.x ~ BETA.y, data=ab)
abline(0,1)
dev.off()

png("../Lifecourse-GWAS-alspac/se_bg_bgen.png")
plot(SE.x ~ SE.y, data=ab)
abline(0,1)
dev.off()


get_lambda <- function(pvector) {
    pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
    o <- -log10(sort(pvector, decreasing=FALSE))
    e <- -log10(ppoints(length(pvector)))
    cs <- qchisq(1-pvector, 1)
    lambda <- median(cs, na.rm=TRUE) / qchisq(0.5, 1)
    return(lambda)
}


get_lambda(10^-ab$LOG10P.x)
get_lambda(10^-ab$LOG10P.y)



lam <- function(fn) {
    res <- readRDS(fn)
    res$pval <- 2 * pnorm(abs(res$BETA / res$SE), lower.tail=FALSE)
    tibble(lambda = get_lambda(res$pval), n=max(res$N), fn=fn)
}

fn <- list.files("/local-scratch/projects/Lifecourse-GWAS/ukb/results2/04", full.names=TRUE) %>% grep(".rds$", ., value=TRUE)


res <- mclapply(fn, lam, mc.cores=50) %>% bind_rows()


vlist <- ab$ID.y

fn1 <- list.files() %>% grep("regenie_test_bestguess", ., value=TRUE)

fn2 <- list.files() %>% grep("regenie_test_bgen", ., value=TRUE)

lam <- function(fn) {
    res <- readRDS(fn)
    res$pval <- 2 * pnorm(abs(res$BETA / res$SE), lower.tail=FALSE)
    tibble(lambda = get_lambda(res$pval), n=max(res$N), fn=fn)
}

fn <- list.files("/local-scratch/projects/Lifecourse-GWAS/ukb/results2/04", full.names=TRUE) %>% grep(".rds$", ., value=TRUE)





c <- readRDS("/local-scratch/projects/Lifecourse-GWAS/ukb/results2/04/bmi_40-45_both.fastGWA.rds")

v <- fread("/local-scratch/projects/Lifecourse-GWAS/ukb/results2/00/variants.txt") %>% as_tibble()

c <- bind_cols(v, c)
c <- subset(c, SNP %in% vlist)
c$pval <- 2 * pnorm(abs(c$BETA / c$SE), lower.tail=FALSE)

get_lambda(c$pval)
