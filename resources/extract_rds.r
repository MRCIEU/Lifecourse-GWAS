library(dplyr)
library(here)

write_phen <- function(rds_file) {
    data <- readRDS(rds_file)
    outdir <- gsub("\\.rds$", "", rds_file)
    dir.create(outdir, showWarnings = FALSE)
    tibble(names(data$age_summary), data$age_summary) %>% write.csv(file.path(outdir, "age_summary.csv"), row.names = FALSE)
    tibble(names(data$age_quantile), data$age_quantile) %>% write.csv(file.path(outdir, "age_quantile.csv"), row.names = FALSE)
    tibble(names(data$deob_summary), data$deob_summary) %>% write.csv(file.path(outdir, "deob_summary.csv"), row.names = FALSE)
    tibble(names(data$sex_table), data$sex_table) %>% write.csv(file.path(outdir, "sex_table.csv"), row.names = FALSE)
    write.csv(data$sex_table, file=file.path(outdir, "sex_table.csv"), row.names = FALSE)
    write.csv(data$sex_yob, file=file.path(outdir, "sex_yob.csv"), row.names = FALSE)
    write.csv(data$categories, file=file.path(outdir, "categories.csv"), row.names = FALSE)
    write.csv(data$sums, file=file.path(outdir, "sums.csv"), row.names = FALSE)
    write.csv(data$categories_ls, file=file.path(outdir, "categories_ls.csv"), row.names = FALSE)
    write.csv(data$sums_ls, file=file.path(outdir, "sums_ls.csv"), row.names = FALSE)
}

write_summary <- function(rds_file) {
    data <- readRDS(rds_file)
    outdir <- gsub("\\.rds$", "", rds_file)
    dir.create(outdir, showWarnings = FALSE)
    write.csv(data$phenotypes, file=file.path(outdir, "phenotypes.csv"), row.names = FALSE)
    tibble(names(data$deob_summary), data$deob_summary) %>% write.csv(file.path(outdir, "deob_summary.csv"), row.names = FALSE)
    tibble(names(data$sex_table), data$sex_table) %>% write.csv(file.path(outdir, "sex_table.csv"), row.names = FALSE)
    write.csv(data$sex_yob, file=file.path(outdir, "sex_yob.csv"), row.names = FALSE)
}

readRenviron(here("config.env"))
resdir <- Sys.getenv("results_dir")


# Get the list of RDS files in 04
files <- list.files(path = file.path(resdir, "04"), pattern = "\\.rds$", full.names = TRUE)

res <- lapply(files, function(file) {
    # Read the RDS file
    readRDS(file)
}) %>% bind_rows()

write.csv(res, file = file.path(resdir, "04-combined_results.csv"), row.names = FALSE)

# Get the lsit of RDS files in 02
files <- list.files(path = file.path(resdir, "02"), pattern = "\\.rds$", full.names = TRUE)
lapply(files, \(x) {
    if(basename(x) == "summary.rds") {
        write_summary(x)
    } else {
        write_phen(x)
    }
})
