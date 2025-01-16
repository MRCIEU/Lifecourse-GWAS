library(here)
library(rmarkdown)

args <- commandArgs(T)
print(args)
rmarkdown::render(args[1], intermediates_dir=args[2], knit_root_dir=args[2], output_dir=args[2])