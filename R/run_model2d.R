# read in arguments
args <- commandArgs(trailingOnly = TRUE)
THE_L <- as.numeric(strsplit(args," ")[[1]][1])

cat("L =",THE_L,"\n")

t1 <- proc.time()
source("R/load_data2.R")
source("R/model_2d.R")
print(proc.time() - t1)

