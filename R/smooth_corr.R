if (!("cdata" %in% ls())) {
	source("R/load_data2.R")
}

# smooth correlations
"smooth_corr" <- function(z, inc=0.005) {
	knots <- seq(min(f), max(f), len=(max(f)-min(f))/inc)

	cor <- vector("list", 9)
	for (i in 1:9) {
		cor[[i]] <- sapply(1:length(knots), function(k) { cor(z[f >= (knots[k]-inc)&f <= (knots[k]+inc),])[i,1:9] })
	}

	list(knots=knots, cor=cor)
}

# smooth covariance
"smooth_cov" <- function(L, z, f, inc=0.02) {
	knots <- seq(min(f), max(f), len=(max(f)-min(f))/inc)
	Nknots <- length(knots)

	Nk <- ncol(z)

	resp.d <- matrix(NA, nrow=Nknots, ncol=Nk)
	resp.o <- matrix(NA, nrow=Nknots, ncol=Nk*(Nk-1)/2)

	for (i in 1:Nknots) {
		Sigma <- cov(z[f >= (knots[i]-inc)&f <= (knots[i]+inc),])

		Sigma.svd  <- svd(Sigma)
		Sigma <- Sigma.svd$u %*% diag(Sigma.svd$d) %*% t(Sigma.svd$u)

		cholSigma <- chol(Sigma)
		resp.d[i,] <- log( diag(cholSigma) )
		resp.o[i,] <- as.vector(cholSigma[upper.tri(cholSigma)])
	}

	Bspline.basis <- create.bspline.basis(c(min(f),max(f)),norder=4,nbasis=L)
	weights <- getbasismatrix(knots, Bspline.basis, nderiv=0)

	ests.d <- matrix(NA, nrow=Nk, ncol=L)
	for (i in 1:Nk) {
		fit <- lm(resp.d[,i] ~ 0+weights)
		ests.d[i,] <- coef(fit)
	}

	ests.o <- matrix(NA, nrow=Nk*(Nk-1)/2, ncol=L)
	for (i in 1:(Nk*(Nk-1)/2)) {
		fit <- lm(resp.o[,i] ~ 0+weights)
		ests.o[i,] <- coef(fit)
	}


	c(as.vector(ests.d), as.vector(ests.o))
}

# weights at specific values of f
"get_weights" <- function(f, L, knots) {
	# construct weights

	if (missing(knots)) {
		#knots <- quantile(f, seq(0,1,length=L))
		knots <- seq(min(f),max(f),length=L)
	}

	if (L > 1) {
		sigma <- 0.5*max(f)/L

		diffs <- sapply(f, function(i) {
			i-knots
		})

		weights <- t(apply(diffs, 1, function(row) {
			exp(-row^2/(2*sigma^2))
		}))

		weights <- apply(weights, 2, function(col) { col/sum(col) })
	} else {
		weights <- matrix(1, nrow=L, ncol=length(f))
	}

	list(w=t(weights), knots=knots)
}

# get some rough least-squares estimates
"ls_estimates_Omega" <- function(sc, L) {
	require(splines)

	if (missing(L)) {
		stop("Must supply L to get LS estimates")
	}

	weights <- bs(f, df=L, intercept=TRUE, Boundary.knots=c(min(f),max(f)))
	knots   <- f[apply(weights, 2, which.max)]
	sc_w    <- predict(weights, sc$knots)

	omega <- rep(NA, L)
	Omega <- vector("list", L)
	for (i in 1:L) {
		Omega[[i]] <- matrix(NA, nrow=9, ncol=9)
	}

	for (i in 1:9) {
		for (j in i:9) {
			fit <- lm(sc$cor[[i]][j,]~0+sc_w)
			for (l in 1:L) {
				Omega[[l]][i,j] <- Omega[[l]][j,i] <- coef(fit)[l]
			}
		}
	}

	# use SVD to then construct p.d. Omega's
	for (l in 1:L) {
		Omega.svd  <- svd(Omega[[l]])
		Omega[[l]] <- cov2cor(Omega.svd$u %*% diag(Omega.svd$d) %*% t(Omega.svd$u))
	}

	Omega
}

# get some rough least-squares estimates
"ls_estimates_cholesky" <- function(data, L) {
	require(fda)

	if (missing(L)) {
		stop("Must supply L to get LS estimates")
	}

	Bspline.basis <- create.bspline.basis(c(min(f),max(f)),norder=4,nbasis=L)
	knots   <- knots(Bspline.basis)
	weights <- getbasismatrix(f, Bspline.basis, nderiv=0)
	sc_w    <- getbasismatrix(sc$knots, Bspline.basis, nderiv=0)
print(head(weights))
print(head(sc_w))
done


	weights <- bs(f, df=L, intercept=TRUE, Boundary.knots=c(min(f),max(f)))
	knots   <- f[apply(weights, 2, which.max)]
	sc_w    <- predict(weights, sc$knots)

	omega <- rep(NA, L)
	Omega <- vector("list", L)
	for (i in 1:L) {
		Omega[[i]] <- matrix(NA, nrow=9, ncol=9)
	}

	for (i in 1:9) {
		for (j in i:9) {
			fit <- lm(sc$cor[[i]][j,]~0+sc_w)
			for (l in 1:L) {
				Omega[[l]][i,j] <- Omega[[l]][j,i] <- coef(fit)[l]
			}
		}
	}

	# use SVD to then construct p.d. Omega's
	for (l in 1:L) {
		Omega.svd  <- svd(Omega[[l]])
		Omega[[l]] <- cov2cor(Omega.svd$u %*% diag(Omega.svd$d) %*% t(Omega.svd$u))
	}

	Omega
}

# smooth correlations
sc <- smooth_corr(z)

if (FALSE) {
	cat("Getting LS estimates for L =",THE_L,"\n")
	# test least squares estimates
	#ests <- ls_estimates_Omega(sc, 5)
	ests <- ls_estimates_cholesky(data, THE_L)
}
