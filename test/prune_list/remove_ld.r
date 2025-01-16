library(data.table)
library(dplyr)

hg19 <- fread("resources/genotypes/hm3_prune_th_hg19.bed.gz") %>% as_tibble()
hg38 <- fread("resources/genotypes/hm3_prune_th_hg38.bed.gz") %>% as_tibble()

a <- subset(hg19, (V1 == 5 & V2 >= 44000000 & V3 <= 51500000) | (V1 == 6 & V2 >= 25000000 & V3 <= 33500000) | (V1 == 8 & V2 >= 8000000 & V3 <= 12000000) | (V1 == 11 & V2 >= 45000000 & V3 <= 57000000))

rsids <- a$V4
hg19 <- subset(hg19, !(V4 %in% rsids))
hg38 <- subset(hg38, !(V4 %in% rsids))

fwrite(hg19, "resources/genotypes/hm3_prune_th_hg19.bed.gz")
fwrite(hg38, "resources/genotypes/hm3_prune_th_hg38.bed.gz")

