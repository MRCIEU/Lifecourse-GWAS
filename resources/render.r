library(here)
library(rmarkdown)

args <- commandArgs(T)
rmarkdown::render(args[1], output_dir=args[2])
