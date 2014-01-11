# get and save initial values for WHICH_CDAT

WHICH_CDAT <- "ST"

source("R/model_spline_cov.R")

for (L in seq(5,20,5)) {
	data <- get_data(L)
	sinits <- smooth_cov(L=L, z=zstar, f=f)
	inits <- get_starts(data, sinits)
	save(sinits, inits, file=paste0("inits/",WHICH_CDAT,"_L",L,".RData"))
}
