library(fda)
library(multicore)
library(MCMCpack)
library(coda)

source("R/load_data2.R")
source("R/smooth_corr.R")
source("R/spline_cov.R")

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
	Bbasis  <- create.bspline.basis(c(min(f),max(f)),norder=4,nbasis=L)
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

	data <- list(prior=100,
		n=n, k=k, y=zstar,
		L=L, weights=weights,
		Nnz=Nnz, Mnz=Mnz-1, Wnz=Wnz,
		uf=uf, ufw=ufw, knots=knots
	)

}

"get_starts" <- function(data, init) {
	# get initial values with BFGS
	if (missing(init)) {
		init <- rep(0, data$L*(data$k+data$k*(data$k-1)/2))
	}

	t1 <- proc.time()
	bfgs <- optim(par=init,
		fn=function(x) {
			lk <- spline_cov_lk(data=data, eval=x)$lk
			-lk
		},
		gr=function(x) {
			gr <- spline_cov_gr(data=data, eval=x)
			-gr
		},
	method="BFGS", control=list(maxit=5000))
	cat("Time to inits: (conv=",bfgs$conv,")\n",sep="")
	print(proc.time()-t1)

	if (bfgs$conv != 0) {
		stop("Bad convergence")
	}

	bfgs$par
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

	if (!has_starts) {
		init <- get_starts(data)
	} else {
		init <- starts
	}

	Nsamples <- round(Niter/thin)

	t1 <- proc.time()
	fits <- mclapply(1:Nchains, mc.cores=Ncores,
		function(i) {
		set.seed(311 + i*Niter);
		fit <- spline_cov(data=data, step_e=step_e, step_L=step_L, thin=thin, inits=init, Niter=Niter, verbose=TRUE)
	})
	cat("Time to samples:\n")
	print(proc.time()-t1)

	# compute DIC

	# get samples, discarding burn
	samples <- vector("list", Nchains)
	sub <- -(1:Nburn)
	res <- vector("list", Nchains)
	dev <- vector("list", Nchains)
	for (i in 1:Nchains) {
		res[[i]] <- matrix(fits[[i]]$samples, nrow=Nsamples)[sub,]
		dev[[i]] <- matrix(fits[[i]]$deviance, nrow=Nsamples)[sub,]
		samples[[i]] <- mcmc( matrix(fits[[i]]$samples, nrow=Nsamples) )
	}
	samples <- mcmc.list(samples)

	# posterior mean
	pmean <- colMeans( do.call(rbind, res) )
	dmean <- mean( unlist(dev) )

	Dbar  <- dmean
	Dtbar <- -2*spline_cov_lk(data=data, eval=pmean)$llik
	pD    <- Dbar - Dtbar
	DIC   <- Dbar + pD

	# save fit
	fname <- paste0("scL",data$L,"_",WHICH_CDAT,".RData")
	save(data, samples, init, DIC, pD, Nburn, thin, Niter, Nsamples, file=paste0("fitsums/fitsum_",fname))

	list(L=data$L, samples=samples, init=init, DIC=DIC, pD=pD, Nburn=Nburn, thin=thin, Niter=Niter, Nsamples=Nsamples)
}

# normal fit
if (exists("WHICH_CDAT") && exists("THE_L")) {
#(c(5,10,15,20)*(9+9*4))^2
#100000/50625
#1.975309*(c(5,10,15,20)*(9+9*4))^2
#	if (WHICH_CDAT == "ST") {
		if (THE_L == 5) {       Niter <- 2*100000;  thin <- Niter/10000; step_e <- 0.10; step_L <- 25; }
		else if (THE_L == 10) { Niter <- 2*400000;  thin <- Niter/10000; step_e <- 0.04; step_L <- 5; }
		else if (THE_L == 15) { Niter <- 2*900000;  thin <- Niter/10000; step_e <- 0.010; step_L <- 5; }
		else if (THE_L == 20) { Niter <- 2*1600000; thin <- Niter/10000; step_e <- 0.00075; step_L <- 10; }
		#else if (THE_L == 15) { step_e <- 0.00001; step_L <- 15; }
		#else if (THE_L == 20) { step_e <- 0.000001; step_L <- 20; }
#	}

	#Niter <- 10000
	#thin <- 10

	Nsamples <- round(Niter/thin)
	Nburn <- round(Nsamples/2)

	print(c(Niter,thin,Nsamples,Nburn))

	#load( paste0("inits/",WHICH_CDAT,"_L",THE_L,".RData") )
	data <- get_data(THE_L)
	inits <- smooth_cov(L=THE_L, z=zstar, f=f)
	fit <- do_fit(data=data, Niter=Niter, Nburn=Nburn, step_e=step_e, step_L=step_L, thin=thin, starts=inits)
	#print(spline_cov_lk(data, inits))
}

if (FALSE) {
Niter <- 20
Nburn <- 5

	#load("inits/ST_L5.RData"); data <- get_data(5)
	#load("inits/ST_L10.RData"); data <- get_data(10)
	#load("inits/ST_L15.RData"); data <- get_data(15)
	#load("inits/ST_L20.RData"); data <- get_data(20)
	#fit5 <- do_fit(data=data, Niter=Niter, Nburn=Nburn, step_e=0.05, step_L=5, starts=inits)
	#fit10 <- do_fit(data=data, Niter=Niter, Nburn=Nburn, step_e=0.01, step_L=5, starts=inits)
	#fit15 <- do_fit(data=data, Niter=Niter, Nburn=Nburn, step_e=0.00001, step_L=10, starts=inits)
	#fit20 <- do_fit(data=data, Niter=Niter, Nburn=Nburn, step_e=0.000001, step_L=1, starts=inits)
}

if (FALSE) {
#fit <- do_fit(0.05, 25) #, good_starts)
#fit <- do_fit(0.025, 5) #, good_starts)
#fit <- do_fit(0.025, 2^9) #, good_starts)

Niter <- 20
Nburn <- 5

data5 <- get_data(5)
init5 <- get_starts(data5)
#fit5 <- do_fit(data=data5, Niter=Niter, Nburn=Nburn, step_e=0.25, step_L=5, starts=init5)

data10 <- get_data(10)
init10 <- get_starts(data10)
#fit10 <- do_fit(data=data10, Niter=Niter, Nburn=Nburn, step_e=0.25, step_L=5, starts=init10)

data15 <- get_data(15)
init15 <- get_starts(data15)
#fit15 <- do_fit(data=data15, Niter=Niter, Nburn=Nburn, step_e=0.25, step_L=5, starts=init15)

data20 <- get_data(20)
init20 <- get_starts(data20)
#fit20 <- do_fit(data=data20, Niter=Niter, Nburn=Nburn, step_e=0.25, step_L=5, starts=init20)

data25 <- get_data(25)
init25 <- get_starts(data25)
#fit25 <- do_fit(data=data25, Niter=Niter, Nburn=Nburn, step_e=0.25, step_L=5, starts=init25)

data30 <- get_data(30)
init30 <- get_starts(data30)
#fit30 <- do_fit(data=data30, Niter=Niter, Nburn=Nburn, step_e=0.25, step_L=5, starts=init30)

#fit1 <- do_fit(get_data(5), 0.025, 5) #, good_starts)
#fit2 <- do_fit(0.025, 10) #, good_starts)
#fit3 <- do_fit(0.025, 15) #, good_starts)
#fit4 <- do_fit(0.025, 20) #, good_starts)

#eps <- .0001; f1 <- do_fit(eps, 10*eps)
#eps <- .00001; f2 <- do_fit(eps, 10*eps)
#eps <- .000001; f3 <- do_fit(eps, 10*eps)
#eps <- .0000001; f4 <- do_fit(eps, 10*eps)
#eps <- .0001; f1 <- do_fit(eps, 10*eps)
#eps <- .0001; f2 <- do_fit(eps, 50*eps)
#eps <- .0001; f3 <- do_fit(eps, 100*eps)
#eps <- .0001; f4 <- do_fit(eps, 150*eps)
#eps <- .0001; f5 <- do_fit(eps, 200*eps)
#eps <- .0001; f6 <- do_fit(eps, 250*eps)
#eps <- .0001; f7 <- do_fit(eps, 512*eps)
#eps <- .0001; f8 <- do_fit(eps, 1024*eps)
#eps <- .0001; fs <- do_fit(eps, 1024*eps)

#round(sapply(1:nrow(ufw),function(i){ 2*invlogit( sum(v3*ufw[i,]) )-1 }),2)
}
