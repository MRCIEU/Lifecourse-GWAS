library(data.table)
library(dplyr)

a <- fread(commandArgs(T)[1])
a$V2 <- make.unique(a$V2)
fwrite(a, commandArgs(T)[1], col.names = F, row.names = F, quote = F, sep = "\t")
