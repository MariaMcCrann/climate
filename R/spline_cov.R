# function to fit the spline covariance model

"spline_cov" <- function(
	prior,
	n, k, y, L, Nnz, Mnz, Wnz,
	step_e, step_L, inits, Niter, samples,
	verbose=FALSE
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
		verbose=as.logical(verbose),
	NAOK=TRUE)

	dyn.unload("RsplineCov.so")

	fit
}

"spline_cov_lk" <- function(
	prior,
	n, k, y, L, Nnz, Mnz, Wnz,
	eval
) {
	dyn.load("RsplineCov.so")

	Mnz <- Mnz-1  # let indices start at 0

	lk <- 0

	res <- .C("spline_cov_lk",
		prior=as.double(prior),
		n=as.integer(n), k=as.integer(k), y=as.double(y),
		L=as.integer(L), Nnz=as.integer(Nnz), Mnz=as.integer(Mnz), Wnz=as.double(Wnz),
		eval=as.double(eval),
		lk=as.double(lk),
	NAOK=TRUE)

	dyn.unload("RsplineCov.so")

	res$lk
}

"spline_cov_gr" <- function(
	prior,
	n, k, y, L, Nnz, Mnz, Wnz,
	eval
) {
	dyn.load("RsplineCov.so")

	Mnz <- Mnz-1  # let indices start at 0

	gr <- rep(0, length(eval))

	res <- .C("spline_cov_gr",
		prior=as.double(prior),
		n=as.integer(n), k=as.integer(k), y=as.double(y),
		L=as.integer(L), Nnz=as.integer(Nnz), Mnz=as.integer(Mnz), Wnz=as.double(Wnz),
		eval=as.double(eval),
		gr=as.double(gr),
	NAOK=TRUE)

	dyn.unload("RsplineCov.so")

	res$gr
}
