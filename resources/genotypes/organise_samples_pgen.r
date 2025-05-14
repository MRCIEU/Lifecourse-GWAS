library(data.table)
library(dplyr)

args <- commandArgs(T)

psam_file <- args[1]
out <- args[2]

all_samples <- fread(psam_file, header=TRUE)
names(all_samples) <- c("V1", "V2", "sex")

message("Unique individuals in sample: ", nrow(all_samples))

if(length(args) == 3) {
    inclusions <- fread(args[3], header=FALSE)
    message("Individuals in inclusion file: ", nrow(inclusions))
    all_samples <- dplyr::filter(all_samples, paste(V1, V2) %in% paste(inclusions$V1, inclusions$V2))
} else {
    message("No inclusion file provided, using all individuals")
}

message("Final number of inclusion individuals: ", nrow(all_samples))

write.table(all_samples, out, row=F, col=F, qu=F, sep=" ")
