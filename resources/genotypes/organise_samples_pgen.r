library(data.table)
library(dplyr)
library(here)
source(here("resources", "genotypes", "read_psam.r"))

args <- commandArgs(T)

psam_file <- args[1]
out <- args[2]

all_samples <- read_psam(psam_file)

message("Unique individuals in sample: ", nrow(all_samples))

if(length(args) == 3) {
    inclusions <- fread(args[3], header=FALSE, keepLeadingZeros=TRUE)
    message("Individuals in inclusion file: ", nrow(inclusions))
    all_samples <- dplyr::filter(all_samples, paste(FAM, IID) %in% paste(inclusions$V1, inclusions$V2))
} else {
    message("No inclusion file provided, using all individuals")
}

message("Final number of inclusion individuals: ", nrow(all_samples))

write.table(all_samples, out, row=F, col=F, qu=F, sep=" ")
