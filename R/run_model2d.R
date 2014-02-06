# read in arguments
args       <- commandArgs(trailingOnly = TRUE)
THE_L      <- as.numeric(strsplit(args," ")[[1]][1])
WHICH_CDAT <- strsplit(args," ")[[1]][2]

cat("L =",THE_L,"\n")
cat("WHICH_CDAT =",WHICH_CDAT,"\n")

t1 <- proc.time()
#source("R/model_2d.R")
#source("R/model_spline_cov.R")
source("R/model_mix_cov.R")
print(proc.time() - t1)

