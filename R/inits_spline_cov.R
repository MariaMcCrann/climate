# get and save initial values for WHICH_CDAT

WHICH_CDAT <- "ST"

source("R/model_spline_cov.R")

for (L in seq(5,30,5)) {
	inits <- get_starts(get_data(L))
	save(inits, file=paste0("inits/",WHICH_CDAT,"_L",L,".RData"))
}
