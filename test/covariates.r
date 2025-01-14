library(dplyr)
library(data.table)

library(here)
library(dotenv)
readRenviron(here("config.env"))

a <- fread(file.path(Sys.getenv("genotype_processed_dir"), paste0("pcs.txt")), header=TRUE)

b <- tibble(FID=a$FID, IID=a$IID, yob = sample(1920:1980, length(FID), replace=TRUE), sex=sample(1:2, length(FID), replace=TRUE), bp_med = rbinom(length(FID), 1, 0.2), choleterol_med= rbinom(length(FID), 1, 0.2), parity=sample(1:5, length(FID), replace=TRUE))