# function to fit the spline covariance model

"spline_cov" <- function(
	data,
	step_e, step_L, inits, Niter,
	verbose=FALSE
) {
	dyn.load("RsplineCov.so")

	fit <- .C("spline_cov_fit",
		prior=as.double(data$prior),
		n=as.integer(data$n), k=as.integer(data$k), y=as.double(data$y),
		L=as.integer(data$L), Nnz=as.integer(data$Nnz), Mnz=as.integer(data$Mnz), Wnz=as.double(data$Wnz),
		step_e=as.double(step_e), step_L=as.integer(step_L),
		inits=as.double(inits),
		Niter=as.integer(Niter),
		samples=as.double(matrix(0, nrow=Niter, ncol=data$L*(data$k + data$k*(data$k-1)/2))),
		deviance=as.double(matrix(0, nrow=Niter, ncol=1)),
		verbose=as.logical(verbose),
	NAOK=TRUE)

	dyn.unload("RsplineCov.so")

	fit
}

"spline_cov_dic" <- function(fit, Nburn=0) {
	if (Nburn > 0) {
		sub <- -(1:Nburn)
	} else {
		sub <- 1:fit$Niter
	}

	# compute DIC
	# posterior mean
	pmean <- colMeans( matrix(fit$samples, nrow=fit$Niter)[sub,] )
	dmean <- mean(fit$deviance[sub])

	Dbar  <- dmean
	Dtbar <- -2*spline_cov_lk(prior=fit$prior, n=fit$n, k=fit$k, y=fit$y, L=fit$L, Nnz=fit$Nnz, Mnz=fit$Mnz, Wnz=fit$Wnz, eval=pmean)$llik
	pD   <- Dbar - Dtbar

	list(DIC=Dbar+pD, Dbar=Dbar, pD=pD)
}

"spline_cov_lk" <- function(data, eval) {
	dyn.load("RsplineCov.so")

	lk   <- 0
	llik <- 0
	lpri <- 0

	res <- .C("spline_cov_lk",
		prior=as.double(data$prior),
		n=as.integer(data$n), k=as.integer(data$k), y=as.double(data$y),
		L=as.integer(data$L), Nnz=as.integer(data$Nnz), Mnz=as.integer(data$Mnz), Wnz=as.double(data$Wnz),
		eval=as.double(eval),
		lk=as.double(lk), llik=as.double(llik), lpri=as.double(lpri),
	NAOK=TRUE)

	dyn.unload("RsplineCov.so")

	list(lk=res$lk, llik=res$llik, lpri=res$lpri)
}

"spline_cov_gr" <- function(data, eval) {
	dyn.load("RsplineCov.so")

	gr <- rep(0, length(eval))

	res <- .C("spline_cov_gr",
		prior=as.double(data$prior),
		n=as.integer(data$n), k=as.integer(data$k), y=as.double(data$y),
		L=as.integer(data$L), Nnz=as.integer(data$Nnz), Mnz=as.integer(data$Mnz), Wnz=as.double(data$Wnz),
		eval=as.double(eval),
		gr=as.double(gr),
	NAOK=TRUE)

	dyn.unload("RsplineCov.so")

	res$gr
}
