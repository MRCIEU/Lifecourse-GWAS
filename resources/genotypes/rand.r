library(data.table)
library(dplyr)

args <- commandArgs(T)

fam <- fread(args[1]) %>% as_tibble()
fam <- dplyr::select(fam, V1, V2) %>%
    mutate(phen=rnorm(n()))

fwrite(fam, args[2], row=F, col=F, qu=F, sep=" ")
