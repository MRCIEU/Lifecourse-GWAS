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

# bmi

traits <- tribble(
    ~trait, ~id,
    "bmi", "ieu-b-40",
    "ldl", "ieu-b-110"
)

lapply(2:nrow(traits), \(i) {
    a <- tophits(traits$id[i]) %>% 
        standardise(., oa_col="nea", pos_col="position") %>%
        select(vid, ea, beta)
    write.table(b, file=here("resources", "genotypes", "tophits", paste0(traits$trait[i], ".txt")), row=F, col=F, qu=F)
})


