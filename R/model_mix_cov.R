library(fda)
library(multicore)
library(MCMCpack)
library(coda)
library(corpcor)

source("R/load_data2.R")
source("R/smooth_corr.R")
#source("R/spline_cov.R")

"fit_mix" <- function(data, prior, Niter=10) {
	# get sqrt() of weights
	rt_w <- sqrt(data$w)

	# get index of max weight
	m_w <- apply(data$w, 1, which.max)

	alpha  <- array(0, dim=c(data$n, data$L, data$k))
	Omega  <- array(NA, dim=c(data$L, data$k, data$k))
	invOmega  <- array(NA, dim=c(data$L, data$k, data$k))

if (TRUE) {
	# initialize alpha to be fraction of observed
	sapply(1:data$n, function(i) {
		sapply(which(data$w[i,] > 0), function(l) {
			if (l != m_w[i]) {
				alpha[i,l,] <<- (data$w[i,l]/(1-data$w[i,m_w[i]]) * data$y[i,])/rt_w[i,l]
			}
		})
	})
}

	keep.Omega  <- array(NA, dim=c(Niter,data$L,data$k,data$k))
	keep.ll     <- rep(NA, Niter)

	# means of each observation
	mu <- matrix(0, nrow=data$n, ncol=data$k)

	for (iter in 1:Niter) {
		# compute mean of each observation
		sapply(1:data$n, function(i) {
			mu[i,] <<- 0
			sapply(which(data$w[i,] > 0), function(l) {
				if (l != m_w[i]) {
					mu[i,] <<- mu[i,] + rt_w[i,l] * alpha[i,l,]
				}
			})
		})

if (TRUE) {
		# update Omega
		sapply(1:data$L, function(l) {
			y_w     <- which(m_w==l)
			pos_w   <- which(data$w[,l] > 0&m_w != l)
			inL     <- length(pos_w) + length(y_w)

			S_y     <- Reduce('+', lapply(y_w,   function(i) { tcrossprod(data$y[i,]-mu[i,])/data$w[i,l] }))
			S_alpha <- Reduce('+', lapply(pos_w, function(i) { tcrossprod(alpha[i,l,]) }))
#print(S_y);print(solve(S_y));
#print(S_alpha);print(solve(S_alpha));
#print(solve(prior$S + S_y + S_alpha))
#print((prior$S + S_y + S_alpha))
#print((prior$S + S_y + S_alpha)/(prior$nu + inL - data$k - 1))
#print(cov2cor( (prior$S + S_y + S_alpha)/(prior$nu + inL - data$k - 1)))
			a       <- prior$nu + inL
			B       <- make.positive.definite(prior$S + S_y + S_alpha)

if (l==1) {
	#print( a-data$k-1 )
	#print(round(S_alpha,3))
	#print(a)
	#print( cov2cor( round(B/(a-data$k-1),3) ) )
}
			Omega[l,,] <<- riwish(a, B)
		})

		# invert Omega's
		sapply(1:data$L, function(l) {
			invOmega[l,,] <<- chol2inv(chol(Omega[l,,]))
		})
}

if (TRUE) {
		# update alpha
		sapply(1:data$L, function(l) {
#F <- (data$w[i,l]/sigma2)*diag(data$k)
#print(chol2inv(chol(F + invOmega)))

			pos_w   <- which(data$w[,l] > 0&m_w != l)

			sapply(pos_w, function(i) {
				d <- rep(0, data$k)
				sapply(which(data$w[i,] > 0), function(k) {
					if (k != m_w[i] && k != l) {
						d <<- d + rt_w[i,k] * alpha[i,k,]
					}
				})

				update_var  <- chol2inv(chol( (data$w[i,l]/data$w[i,m_w[i]])*invOmega[m_w[i],,] + invOmega[l,,]))
				update_mu   <- (rt_w[i,l]/data$w[i,m_w[i]])*(update_var %*% invOmega[m_w[i],,] %*% (data$y[i,]-d) )
				alpha[i,l,] <<- update_mu + t(chol(update_var)) %*% rnorm(data$k)
			})
		})
}

		ll <- sum(sapply(1:data$n, function(i) {
			S <- Reduce('+', lapply(which(data$w[i,] > 0), function(l) { data$w[i,l]*Omega[l,,] }))
			cholS <- chol(S)
			invS <- chol2inv(cholS)
			-sum(log(diag(cholS))) -0.5*t(data$y[i,]) %*% invS %*% (data$y[i,])
		}))

		if (iter %% 10==0) {
			print(iter)
			print(ll)
			print( round( cov2cor( Omega[1,,] ), 3) )
		}

		keep.Omega[iter,,,] <- Omega
		keep.ll[iter]       <- ll
	}

	list(Omega=keep.Omega, ll=keep.ll)
}

if (FALSE) {
	# let's subset the data...
	z <- z[,c(1,5,9)]

	#keep <- T <= 40
	#keep <- c(1, 1+sort( sample.int(nrow(z)-1, size=round(nrow(z)/8)) ))
	#keep <- c(1, round(seq(2, nrow(z), len=round(nrow(z)/8))) )
	#z <- z[keep,]; d <- d[keep]; f <- f[keep]; T <- T[keep]
}

# normalize the data
zstar <- sqrt(d) * z

# get range of f before we remove d_i=0
range.f <- range(f)

# remove d_i = 0
rem <- which(d==0)
if (length(rem) > 0) {
	z <- z[-rem,]; d <- d[-rem]; f <- f[-rem]; T <- T[rem]; zstar <- zstar[-rem,]
}

n  <- nrow(zstar)
k  <- ncol(zstar)

"get_data" <- function(L) {
	if (L < 4) { stop("B-splines require L > 3\n") }

	# create basis functions
	#Bbasis  <- create.bspline.basis(c(min(f),max(f)),norder=4,nbasis=L)
	Bbasis  <- create.bspline.basis(range.f,norder=4,nbasis=L)
	knots   <- knots(Bbasis)
	weights <- getbasismatrix(f, Bbasis, nderiv=0)
	uf      <- quantile(f, seq(0,1,length=100))
	ufw     <- getbasismatrix(uf, Bbasis, nderiv=0)

	# capture which are non-zero
	nz <- apply(weights, 1, function(x){ which(x!=0) })

	# get number of non-zeros
	Nnz <- sapply(1:length(nz), function(i){ length(nz[[i]]) })

	# get non-zero indices
	Mnz <- matrix(0, nrow=length(nz), ncol=max(Nnz))
	sapply(1:length(nz), function(i){ Mnz[i,1:Nnz[i]] <<- nz[[i]] })

	# get non-zero weights
	Wnz <- matrix(0, nrow=length(nz), ncol=max(Nnz))
	sapply(1:length(nz), function(i){ Wnz[i,1:Nnz[i]] <<- weights[i,nz[[i]]] })

#set.seed(311)
#fit <- fit_mix(data=list(y=y, n=n, k=k, L=L, w=weights),
#	prior=list(nu=10, S=(10-k-1)*diag(k), alpha_s=0.001, beta_s=0.001)
#)

	data <- list(prior=list(nu=20, S=(20-k-1)*diag(k)),
		n=n, k=k, y=zstar, f=f,
		L=L, w=weights,
		Nnz=Nnz, Mnz=Mnz-1, Wnz=Wnz,
		uf=uf, ufw=ufw, knots=knots,
		sigma2=0.01
	)

}

"do_fit" <- function(data, Niter=100, Nburn=50, step_e=0.01, step_L=1, thin=1, starts) {
	# do we have starting values?
	if (missing(starts))
		has_starts <- FALSE
	else
		has_starts <- TRUE

	Nchains <- 3
	Ncores  <- 3
	Nparam <- data$L*(data$k + data$k*(data$k-1)/2)

if (FALSE) {
	if (!has_starts) {
		init <- get_starts(data)
	} else {
		init <- starts
	}
}

	Nsamples <- round(Niter/thin)

	t1 <- proc.time()
	#fits <- mclapply(1:Nchains, mc.cores=Ncores,
	fits <- lapply(1:Nchains,
		function(i) {
		set.seed(311 + i*Niter);
		fit <- fit_mix(data=data, prior=list(nu=10, S=(10-k-1)*diag(k)), Niter=Niter)
	})
	cat("Time to samples:\n")
	print(proc.time()-t1)

if (TRUE) {
	# compute DIC

	# get samples, discarding burn
	sub <- -(1:Nburn)
	Omegas <- vector("list", Nchains)
	dev <- vector("list", Nchains)
	for (i in 1:Nchains) {
		Omegas[[i]] <- fits[[i]]$Omega[sub,,,]
		dev[[i]] <- -2*fits[[i]]$ll[sub]
	}

	# posterior mean
	mOmegas <- array(NA, dim=c(data$L, data$k, data$k))
	sapply(1:data$L, function(l) {
		sapply(1:data$k, function(k1) {
			sapply(1:data$k, function(k2) {
				mOmegas[l,k1,k2] <<- mean(c(Omegas[[1]][,l,k1,k2],Omegas[[2]][,l,k1,k2],Omegas[[3]][,l,k1,k2]))
			})
		})
	})
print(mOmegas[1,,])
print(mOmegas[2,,])

	Dbar  <- mean(c(dev[[1]],dev[[2]],dev[[3]]))
	Dtbar <- -2*sum(sapply(1:data$n, function(i) {
			S <- Reduce('+', lapply(which(data$w[i,] > 0), function(l) { data$w[i,l]*mOmegas[l,,] }))
			cholS <- chol(S)
			invS <- chol2inv(cholS)
			-sum(log(diag(cholS))) -0.5*t(data$y[i,]) %*% invS %*% (data$y[i,])
		}))
	pD    <- Dbar - Dtbar
	DIC   <- Dbar + pD

	# save fit
	fname <- paste0("mcL",data$L,"_",WHICH_CDAT,".RData")
	#save(data, samples, init, DIC, pD, Nburn, thin, Niter, Nsamples, file=paste0("fitsums/fitsum_",fname))
	save(data, mOmegas, Omegas, DIC, pD, Nburn, thin, Niter, file=paste0("fitsums/fitsum_",fname))
}

	#list(L=data$L, samples=samples, init=init, DIC=DIC, pD=pD, Nburn=Nburn, thin=thin, Niter=Niter, Nsamples=Nsamples)
	list(fits=fits, L=data$L, Omegas=mOmegas, DIC=DIC, pD=pD, Nburn=Nburn, thin=thin, Niter=Niter)
	#list(fits)
}

# normal fit
if (exists("WHICH_CDAT") && exists("THE_L")) {
	Niter <- 1000
	thin <- 1

	Nsamples <- round(Niter/thin)
	Nburn <- round(Nsamples/2)

	print(c(Niter,thin,Nsamples,Nburn))

	cat("Getting data\n")
	data <- get_data(THE_L)

	cat("Getting inits\n")
	#inits <- smooth_cov(L=THE_L, z=data$y, f=data$f)

	cat("Running fit\n")
	fit <- do_fit(data=data, Niter=Niter, Nburn=Nburn)
	#print(spline_cov_lk(data, inits))
}
