library(dplyr)
library(data.table)


pl <- fread("phenotype_list.csv") %>% as_tibble()

pl %>% as.data.frame

covs <- lapply(1:nrow(pl), \(i) {
    tibble(
        trait=pl$pheno_id[i],
        covs=strsplit(pl$covs[i], split=":") %>% unlist
    )
}) %>% bind_rows() %>% filter(covs != "sex")

ids <- sample(1000:9999, 20)
age <- (sample((20*12):(80*12), 20) / 12) %>% round(., 2)


for(i in 1:nrow(pl)) {
    value <- runif(length(ids), min=pl$min[i], max=pl$max[i]) %>% round(., 2)
    dat <- tibble(FID=ids, IID=ids, value, age)
    dat2 <- dat[1:10, ]
    dat2$age <- dat2$age + 5
    dat2$value <- runif(nrow(dat2), min=pl$min[i], max=pl$max[i]) %>% round(., 2)
    dat <- bind_rows(dat2, dat)
    co <- subset(covs, trait==pl$pheno_id[i])$covs
    if(length(co) > 0) {
        dat <- dat %>% mutate(!!co := sample(0:1, nrow(dat), replace=TRUE))
    }
    dat <- arrange(dat, FID)
    write.table(dat, file=file.path("test", "example-data", "phen_input", paste0(pl$pheno_id[i], ".txt")), sep="\t", quote=FALSE, row.names=FALSE)
}

dat


sc <- tibble(FID=ids, IID=ids, sex=sample(1:2, length(ids), replace=TRUE), yob=sample(1920:1970, length(ids), replace=TRUE), batch=sample(1:2, length(ids), replace=TRUE))
write.table(sc, file=file.path("test", "example-data", "phen_input", "static_covariates.txt"), sep="\t", quote=FALSE, row.names=FALSE)


