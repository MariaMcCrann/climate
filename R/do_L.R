# run model for a specific L
THE_L <- as.numeric(read.table("L")[1])

cat("Running model for",THE_L,"\n")
source("R/load_data2.R")
source("R/model_2d.R")
