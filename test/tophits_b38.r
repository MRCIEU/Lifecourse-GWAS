library(dplyr)
library(glue)
library(data.table)
library(tidyr)
library(here)


update_tophit_build <- function(tophit_file, output_file, liftover_binary=here("test/prune_list/liftOver"), chain_file=here("test/prune_list/hg19ToHg38.over.chain")) {
    a <- fread(tophit_file) %>% as_tibble() %>% separate(V1, remove=FALSE, into=c("chr", "pos", "ea", "oa"))
    bed <- tibble(chr=paste0("chr",a$chr), p1=a$pos, p2=a$pos, id=a$V1) %>% mutate(chr=gsub("chr23", "chrX", chr))

    tf <- tempfile()
    write.table(bed, file=tf, row=F, col=F, qu=F)

    tf

    cmd <- glue("{liftover_binary} {tf} {chain_file} {tf}.lo {tf}.unmapped")
    system(cmd)

    b <- fread(glue("{tf}.lo"), header=FALSE) %>% as_tibble()
    b <- left_join(a, b, by=c("V1"="V4")) %>% mutate(newid=paste0(chr, ":", V2.y, "_", ea, "_", oa))
    b$newid

    a
    s <- tibble(b$newid, b$V2.x, b$V3.x)
    a

    write.table(s, file=output_file, row=F, col=F, qu=F)
}

tophit_file <- "/home/gh13047/repo/Lifecourse-GWAS-ukb/resources/genotypes/tophits/adiponectin.txt"

dir.create(here("resources/genotypes/tophits/hg38"))

list.files(here("resources/genotypes/tophits/hg19")) %>% lapply(., \(x) {
    print(x)
    update_tophit_build(here("resources/genotypes/tophits/hg19", x), here("resources/genotypes/tophits/hg38", x))
})

