library(data.table)
library(dplyr)

args <- commandArgs(T)

phenolist_file <- args[1]
# phenolist_file <- "/local-scratch/projects/Lifecourse-GWAS/ukb/phen_proc2/phenolist"

dir.create(file.path(dirname(phenolist_file), "regenie"), showWarnings = F)

phenolist <- scan(phenolist_file, what="character") %>% grep("bmi", ., value=T)

dat <- tibble(fn=phenolist, phen=basename(fn) %>% gsub(".phen", "", .))
a <- fread(phenolist[1], he=F)
names(a) <- c("FID", "IID", dat$phen[1])
for(i in 2:nrow(dat)) {
  message(i)
  b <- fread(phenolist[i], he=F)
  names(b) <- c("FID", "IID", dat$phen[i])
  a <- full_join(a, b, by=c("FID", "IID"))
}
phen <- a

dat2 <- dat
dat2$fn <- gsub(".phen$", ".covs", dat2$fn)

covs <- lapply(1:nrow(dat2), \(i) {
  message(i)
  b <- fread(dat2$fn[i], he=F)
  names(b)[1:2] <- c("FID", "IID")
  b
}) %>% bind_rows() %>% filter(!duplicated(paste(FID, IID)))

write.table(phen, file = file.path(dirname(phenolist_file), "regenie", "phen.txt"), quote = F, row.names = F, sep=" ", col.names = T)
write.table(covs, file = file.path(dirname(phenolist_file), "regenie", "covs.txt"), quote = F, row.names = F, sep=" ", col.names = T)

# Static covariates that are the same for all phenotypes are in the covs file
# Other covariates (e.g. time-varying ones) need to be adjusted from the phenotype at the phenotype organisation stage.

