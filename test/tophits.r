library(ieugwasr)
library(dplyr)
library(here)

standardise <- function(d, ea_col="ea", oa_col="oa", beta_col="beta", eaf_col="eaf", chr_col="chr", pos_col="pos", vid_col="vid") {
    toflip <- d[[ea_col]] > d[[oa_col]]
    d[[eaf_col]][toflip] <- 1 - d[[eaf_col]][toflip]
    d[[beta_col]][toflip] <- d[[beta_col]][toflip] * -1
    temp <- d[[oa_col]][toflip]
    d[[oa_col]][toflip] <- d[[ea_col]][toflip]
    d[[ea_col]][toflip] <- temp
    d[[vid_col]] <- paste0(d[[chr_col]], ":", d[[pos_col]], "_", d[[ea_col]], "_", d[[oa_col]])
    d
}

traits <- read.csv(here("phenotype_list.csv"))
traits <- bind_rows(traits, tibble(Name = "Pulse pressure", "pheno_id" = "pp", "opengwasid" = "ebi-a-GCST90018970"))
i <- 23
lapply(1:nrow(traits), \(i) {
    of <- here("resources", "genotypes", "hg19", "tophits", paste0(traits$pheno_id[i], ".txt"))
    if(traits$pheno_id[i] == "depression") {
        return(NULL)
    }
    if(!file.exists(of)) {
        message(i, " " , of)
        a <- try(tophits(traits$opengwasid[i]))
        if(! "try-error" %in% class(a)) 
        {
            if(nrow(a) > 0) {
                a <- a %>% 
                    standardise(., oa_col="nea", pos_col="position") %>%
                    select(vid, ea, beta)
                write.table(a, file=of, row=F, col=F, qu=F)
            }
            if(traits$pheno_id[i] == "bmi") {
                of <- here("resources", "genotypes", "hg19", "tophits", paste0("bmiz.txt"))
                write.table(a, file=of, row=F, col=F, qu=F)
            }
        }
    }
})



