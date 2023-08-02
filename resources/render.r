library(here)
library(rmarkdown)
library(jsonlite)

args <- commandArgs(T)

config <- jsonlite::read_json(here("config.json"))
rmarkdown::render(args[1], output_dir=config$results_dir)
