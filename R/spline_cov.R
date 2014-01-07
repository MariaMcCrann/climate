# function to fit the spline covariance model

"spline_cov" <- function(
	prior,
	n, k, y, L, Nnz, Mnz, Wnz,
	step_e, step_L, inits, Niter, samples
) {
	dyn.load("RsplineCov.so")

	Mnz <- Mnz-1  # let indices start at 0

	fit <- .C("spline_cov_fit",
		prior=as.double(prior),
		n=as.integer(n), k=as.integer(k), y=as.double(y),
		L=as.integer(L), Nnz=as.integer(Nnz), Mnz=as.integer(Mnz), Wnz=as.double(Wnz),
		step_e=as.double(step_e), step_L=as.integer(step_L),
		inits=as.double(inits),
		Niter=as.integer(Niter), samples=as.double(samples),
	NAOK=TRUE)

	dyn.unload("RsplineCov.so")

	fit
}
