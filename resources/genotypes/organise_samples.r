library(data.table)
library(dplyr)

args <- commandArgs(T)

geno_files <- args[1]
out <- args[2]

file_list <- fread(geno_files) %>% as_tibble()
all_samples <- lapply(1:nrow(file_list), \(i) {
    fread(file_list$V2[i], skip=2, header=FALSE) %>% as_tibble() %>%
        dplyr::select(V1, V2)
}) %>% 
    bind_rows() %>%
    dplyr::filter(!duplicated(paste(V1, V2)))

message("Unique individuals in sample: ", nrow(all_samples))

if(length(args) == 3) {
    inclusions <- fread(args[3], header=FALSE)
    message("Individuals in inclusion file: ", nrow(inclusions))
    all_samples <- dplyr::filter(all_samples, paste(V1, V2) %in% paste(inclusions$V1, inclusions$V2))
}

message("Final number of inclusion individuals: ", nrow(all_samples))

fwrite(all_samples, out, row=F, col=F, qu=F, sep=" ")
