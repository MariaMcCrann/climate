require(fields)
require(maps)
library(multicore)
library(rstan); set_cppo("fast")
library(MCMCpack)
library(mvtnorm)

source("R/plot_2d.R")
source("R/load_data2.R")

if (TRUE) { # work with data

# plot this data
if (TRUE) {pdf("pdf/2d/data.pdf")
	par(bty="l")
	par(mfrow=c(2,2))
	image.plot(tlon,tlat,temps[,,9],xlab="",ylab="",axes=F,main="Data");   map("world",add=T)
	image.plot(tlon,tlat,temps[,,1],xlab="",ylab="",axes=F,main="GCM");    map("world",add=T)
	image.plot(tlon,tlat,temps[,,5],xlab="",ylab="",axes=F,main="RCM"); map("world",add=T)
	#image.plot(tlon,tlat,temps[,,6],xlab="",ylab="",axes=F,main="RCM(6)"); map("world",add=T)
graphics.off()}

# plot Zs
if (FALSE) {pdf("pdf/2d/data_trans.pdf")
	par(bty="l")
	par(mfrow=c(2,2))
  image.plot(tlon,tlat,Z[,,9],xlab="",ylab="",axes=F,main="Data")
  image.plot(tlon,tlat,Z[,,1],xlab="",ylab="",axes=F,main="GCM")
  image.plot(tlon,tlat,Z[,,2],xlab="",ylab="",axes=F,main="RCM(2)")
  image.plot(tlon,tlat,Z[,,6],xlab="",ylab="",axes=F,main="RCM(6)")
graphics.off()}

# let's fit a model
dim2_code <- "
	data {
		int<lower=0>  n;         // sample size
		int<lower=0>  k;         // number of data sources
		matrix[n,k]   Z;         // data Z
		vector[n]     d;         // diag(D) from decomposition of Q

		int<lower=0>  Nt;        // how many T_i do we have?
		int<lower=0>  T[n];      // time T_i for each Z

		int<lower=0>  L;         // number of cov matrices in mixture
		matrix[L,Nt]   weights;   // resolution weights for each T
	}
	transformed data {
		matrix[n,k] Zstar;
		vector[k]   zero[n];

		// normalize Z
		for (i in 1:n)
			Zstar[i] <- sqrt(d[i])*Z[i];

		// mean for the delta's
		for (i in 1:n)
				for (j in 1:k)
					zero[i,j] <- 0.0;
	}
	parameters {
		vector<lower=0>[k]           omega[L];
		corr_matrix[k]               corrOmega[L];
	}
	model {
		matrix[k,k] Omega[L];
		matrix[k,k] Sigma_d[Nt];

		// construct Omega
		for (l in 1:L)
			for (i in 1:k)
				for (j in 1:k)
					Omega[l,i,j] <- omega[l,i]*omega[l,j]*corrOmega[l,i,j];

		// construct Sigma(T) = sum_l w_l(T) Omega_l
		for (t in 1:Nt) {
			Sigma_d[t] <- weights[1,t] * Omega[1];

			for (l in 2:L)
				Sigma_d[t] <- Sigma_d[t] + weights[l,t] * Omega[l];
		}

		// priors
		for (l in 1:L)
			omega[l] ~ cauchy(0, 5);

		// model
		for (i in 2:n)
			Zstar[i]' ~ multi_normal(zero[i], Sigma_d[T[i]]);
	}

	generated quantities {
#		matrix[k,k] Omega[L];
#		matrix[k,k] Sigma_d[Nt];
		matrix[k,k] corrSigma_d[Nt];

#		// construct Omega
#		for (l in 1:L)
#			for (i in 1:k)
#				for (j in 1:k)
#					Omega[l,i,j] <- omega[l,i]*omega[l,j]*corrOmega[l,i,j];
#
#		// construct Sigma(T) = sum_l w_l(T) Omega_l
#		for (t in 1:Nt) {
#			Sigma_d[t] <- weights[1,t] * Omega[1];
#
#			for (l in 2:L)
#				Sigma_d[t] <- Sigma_d[t] + weights[l,t] * Omega[l];
#		}
#
#		// construct corrSigma(T)
#		for (t in 1:Nt)
#			for (i in 1:k)
#				for (j in 1:k)
#						corrSigma_d[t,i,j] <- Sigma_d[t,i,j]/sqrt(Sigma_d[t,i,i]*Sigma_d[t,j,j]);

		// construct Sigma(T) = sum_l w_l(T) Omega_l
		for (t in 1:Nt) {
			corrSigma_d[t] <- weights[1,t] * corrOmega[1];

			for (l in 2:L)
				corrSigma_d[t] <- corrSigma_d[t] + weights[l,t] * corrOmega[l];
		}

	}
"

"rInvWish" <- function(v, invS) {
	chol2inv(chol(rwish(v, invS)))
}

"model_2d_gibbs" <- cmpfun( function(dat, Niter=1000) {
	n <- dat$n
	k <- dat$k
	w <- dat$weights
	w2 <- w^2
	L <- dat$L
	T <- dat$T

	# transform data
	Zstar <- sqrt(d) * dat$Z;

	# setup priors
	nu <- k+2
	S  <- diag(k)
	update_nu <- nu+n

	# where do we store current values?
	Omega <- array(0, dim=c(L,k,k))
	cholOmega <- array(0, dim=c(L,k,k))
	if (L > 1) { alpha <- array(0, dim=c(n,L-1,k)) }

	# values at each iter
	i_Omega <- array(0, dim=c(Niter,L,k,k))
	if (L > 1) { i_alpha <- array(0, dim=c(Niter,n,L-1,k)) }

	# get initial values
	for (l in 1:L) {
		Omega[l,,]     <- riwish(nu, S)
		cholOmega[l,,] <- t(chol(Omega[l,,]))
	}

	if (L > 1) {
		alpha <- array(rnorm(n*(L-1)*k), dim=c(n,L-1,k))
	}

	mu <- matrix(0, nrow=n, ncol=k)
	if (L == 1) {
		S_mu <- matrix(0, nrow=k, ncol=k)
 		sapply(1:n, function(i) {
			S_mu <<- S_mu + (Zstar[i,] %*% t(Zstar[i,])) /w[1,T[i]]^2
		})

		update_S <- chol2inv(chol(S + S_mu))
	}

	for (iter in 1:Niter) {
		# update Omega_1
		if (L > 1) {
			for (i in 1:n) {
				for (l in 2:L) {
					mu[i,] <- mu[i,] + w[l,T[i]] * (cholOmega[l,,] %*% alpha[i,l-1,])
				}
			}

			S_mu <- matrix(0, nrow=k, ncol=k)
	 		sapply(1:n, function(i) {
				S_mu <<- S_mu + ( (Zstar[i,]-mu[i]) %*% t(Zstar[i,]-mu[i]) )/w[1,T[i]]^2
			})

			update_S <- chol2inv(chol(S + S_mu))
		}

		Omega[1,,] <- rInvWish(update_nu, update_S)
		cholOmega[1,,] <- t(chol(Omega[1,,]))

		if (L > 1) {
			# update Omega_2, ..., Omega_L
			for (l in 2:L) {
				S_theta <- matrix(0, nrow=k, ncol=k)
		 		sapply(1:n, function(i) {
					theta <- cholOmega[l,,] %*% alpha[i,l-1,]
					S_theta <<- S_theta + tcrossprod(theta)
				})

				Omega[l,,]     <- rInvWish(update_nu, chol2inv(chol(S + S_theta)))
				cholOmega[l,,] <- t(chol(Omega[l,,]))
			}

			# update alpha
			invOmega1 <- chol2inv(t(cholOmega[1,,]))
			M1 <- cholOmega[1,,] %*% invOmega1
			M2 <- M1 %*% t(cholOmega[1,,])
			M3 <- diag(k)
			talpha <- array(rnorm(n*(L-1)*k), dim=c(n,L-1,k))
			for (i in 1:n) {
				for (l in 1:(L-1)) {
					invA <- chol2inv(chol(w2[l+1,T[i]] * M2 + M3))
					b <- w[l+1,T[i]] * M1 %*% Zstar[i,]
					alpha[i,l,] <- invA %*% b + chol(invA) %*% talpha[i,l,]
				}
			}

			i_alpha[iter,,,] <- alpha
		}

		i_Omega[iter,,,] <- Omega
	}

	# compute corr mat Gamma at each time point
	uT <- sort(unique(T))
	nT <- length(uT)

	Gamma <- array(0, dim=c(Niter,nT,k,k))
	for (iter in 1:Niter) {
		for (t in uT) {
			Sigma <- matrix(0, nrow=k, ncol=k)
			for (l in 1:L) {
				Sigma <- Sigma + w[l,t] * i_Omega[iter,l,,]
			}
			Gamma[iter,t,,] <- cov2cor(Sigma)
		}
	}

	if (L > 1) {
		list(Gamma=Gamma, Omega=i_Omega, alpha=i_alpha)
	} else {
		list(Gamma=Gamma, Omega=i_Omega)
	}
})

"get_weights" <- function(d, L) {
	# construct weights
	if (L > 1) {
		knots <- quantile(d, seq(0,1,length=L))
		sigma <- 0.5*max(d)/L

		diffs <- sapply(d, function(i) {
			i-knots
		})

		weights <- t(apply(diffs, 1, function(row) {
			exp(-row^2/(2*sigma^2))
		}))

		weights <- apply(weights, 2, function(col) { col/sum(col) })
	} else {
		weights <- matrix(1, nrow=L, ncol=length(d))
	}

	t(weights)
}

if (TRUE) {
# let's subset the data...
z <- z[,c(1,5,9)]
#keep <- T <= 40
keep <- c(1, 1+sort( sample.int(nrow(z)-1, size=round(nrow(z)/4)) ))
z <- z[keep,]; d <- d[keep]; T <- T[keep]
}

zstar <- sqrt(d) * z

#n  <- n1*n2
n  <- nrow(z)
k  <- ncol(z)
Nt <- max(T)

# compile once...
weights <- get_weights(d, 1); max_w <- rep(1, length(d)); i_max_w <- max_w; ud <- quantile(d, seq(0,1,length=10)); udw <- get_weights(ud, 1)
dat = list(n=n, k=k, Zstar=zstar, L=1, weights=weights, i_max_w=i_max_w, Nud=length(ud), udw=udw)
#fit2d <- stan(file = 'stan/model_2d_re.stan', data = dat, iter = 10, chains = 1,pars="Omega");save(fit2d, file="fit2d.RData");done
#fit2d <- stan(file = 'stan/model_2d.stan', data = dat, iter = 10, chains = 1,pars="Omega");save(fit2d, file="fit2d.RData");done
load("fit2d.RData")

#Ls <- 10 #1:6
Ls <- THE_L
fits <- lapply(Ls, function(L) {
	weights <- get_weights(d, L)
	max_w   <- apply(weights, 1, max)
	i_max_w <- apply(weights, 1, which.max)
	ud      <- quantile(d, seq(0,1,length=100))
	udw     <- get_weights(ud, L)

	dat = list(n=n, k=k, Zstar=zstar,
		L=L, weights=weights, i_max_w=i_max_w,
		Nud=length(ud), udw=udw
	)

fn.inits <- function() {
	system.time(311)

	dq <- quantile(d, probs=(0:L)/L)
	omega <- array(NA, dim=c(L,k))
	corrOmega <- array(NA, dim=c(L,k,k))
	for (l in 1:L) {
		#corrOmega[l,,] <- rWishart(1, k+10, diag(k))[,,1]
		#omega[l,] <- 2*rep(1, k)/l
		omega[l,] <- apply(zstar[d > dq[l]&d <= dq[l+1],], 2, sd)
		corrOmega[l,,] <- diag(k)
		#corrOmega[l,,] <- cov2cor(corrOmega[l,,])
	}

	list("omega"=omega, "corrOmega"=corrOmega, "alpha"=array(rnorm(n*L*k), dim=c(n,L,k)))
}

if (FALSE) { # use gibbs
	fit <- model_2d_gibbs(dat)
} else { # use stan
	delta  <- 0.15; max_td <- 5

if (FALSE) {
	# get a starting point
	cat("Getting starting inits...\n")
	inits <- list(fn.inits())

	cat("Getting starts...\n")
	sflist <- mclapply(1:3, mc.cores=3, function(i) {
		stan(fit=fit2d, data=dat, iter=200, chains=1, chain_id=i, seed=03101983,
		     delta=delta, max_treedepth=max_td, refresh=1, verbose=TRUE, init=inits, pars=c("omega","corrOmega"))
	})
	cat("Merging starts...\n")
	fit.start <- sflist2stanfit(sflist)
	cat("Getting start summary...\n")
	fs.start <- summary(fit.start, probs=c(0.5) )$summary

"fn.start_ini" <- function() {
	list(
		"omega"=aperm(array(fs.start[grep("^omega",rownames(fs.start)),4],dim=c(k,L)), c(2,1)),
		"corrOmega"=aperm( array(fs.start[grep("^corrOmega",rownames(fs.start)),4],dim=c(k,k,L)), c(3,2,1) ),
		"alpha"=inits[[1]]$alpha
#		"alpha"=aperm( array(fs.start[grep("^alpha",rownames(fs.start)),4],dim=c(k,L,n)), c(3,2,1) )
	)
}

	cat("Getting start inits...\n")
	inits <- list( fn.start_ini() )
print(inits[[1]]$omega)
print(inits[[1]]$corrOmega[1,,])
print(inits[[1]]$alpha[1,,])

	print(object.size(fit.start),units="Mb")
	print(object.size(fs.start),units="Mb")
	rm("fit.start")
	rm("fs.start")
}

	# run in parallel
	delta  <- 0.15; max_td <- 6
	Nfits <- 3
	sflist <- mclapply(1:Nfits, mc.cores=3,
	#sflist <- lapply(1:3,
		function(i) {
			tf <- stan(fit=fit2d, data=dat, iter=1000, #init=inits, #epsilon=3.14821e-06,
			     delta=delta, max_treedepth=max_td,
				   #max_tree_depth=6, #nondiag_mass=TRUE, #,delta=0.25,max_tree_depth=-1,equal_step_sizes=TRUE, epsilon=
			     chains = 1, seed=03101983, chain_id=i, refresh=1, verbose=FALSE,
			     #pars=c("sigma_e","omega","corrOmega","corrSigma_d"))
			     #pars=c("omega","corrSigma_d")
			     #pars=c("omega","corrOmega")
			     pars=c("corrSigma_d","Omega")
			)

#			cat("Saving",i,"\n")
#			save(tf, file=paste0("~/tmp/fits/climate",i,".RData"))
#
#			NA
			tf
		})

	if (length(sflist) > 1) {
#		for (i in 1:Nfits) {
#			cat("Loading",i,"\n")
#			load(file=paste0("~/tmp/fits/climate",i,".RData"))
#			sflist[[i]] <- tf
#		}

		cat("Merging fits...\n")
		fit <- sflist2stanfit(sflist)
	} else {
#		load(file=paste0("~/tmp/fits/climate",1,".RData"))
#		fit <- tf
		fit <- sflist[[1]]
	}
	print(object.size(sflist),units="Mb")
	rm("sflist")

	cat("Computing fit summary..\n")
	fitsum <- summary(fit)$summary

	# compute DIC, LPML, GG, and CRPS
	cat("Computing DIC, LPML, GG, and CRPS...\n")
	svars <- list(
		"Omega"=aperm( array(fitsum[grep("^Omega",rownames(fitsum)),1],dim=c(k,k,L)), c(3,2,1) )
	)

print(dat$Zstar[2,])
print(weights[2,])
print(svars$Omega[1,,])

if (FALSE) {
	Omegas <- extract(fit)$Omega
	incr  <- round(nrow(Omegas)/Nfits)

	cat("Computing log lik at each iteration...\n")
	ll_calcs <- unlist( mclapply(1:Nfits, function(block) {
			start <- 1 + (block-1)*incr
			end   <- start + incr
			if (end > nrow(Omegas)) { end <- nrow(Omegas) }

			sapply(start:end, function(iter) {
				ll <- sum( sapply(2:n, function(i) {
					i_Omega <- matrix(0, nrow=k, ncol=k)
					for (l in 1:L) {
						i_Omega <- i_Omega + weights[i,l] * Omegas[iter,l,,]
					}
					dmvnorm(dat$Zstar[i,], rep(0, k), i_Omega, log=TRUE)
				}) )

				ll
			})
		}) )

	cat("Computing DIC...\n")
	dbar <- -2 * mean(ll_cals)
	dhat <- sum( sapply(2:n, function(i) {
		i_Omega <- matrix(0, nrow=k, ncol=k)
		for (l in 1:L) {
			i_Omega <- i_Omega + weights[i,l] * svars$Omega[l,,]
		}
		-2 * dmvnorm(dat$Zstar[i,], rep(0, k), i_Omega, log=TRUE)
	}) )

	pD  <- dbar - dhat
	DIC <- dbar + pD

	cat("Computing LPML...\n")
	LPML <- sum(-log( sum( (1/ll_calcs)/nrow(Omegas) ) ))

	cat("Computing GG...\n")
	GG <- 0
	cat("Computing CRPS...\n")
	CRPS <- 0
} else {
	ll_calcs <- 0
	DIC <- 0
	LPML <- 0
	GG <- 0
	CRPS <- 0
}

	print(object.size(fit),units="Mb")
	rm("fit")

#	# make plots
#	try({
#		print(summary(fitsum$summary[,10]))
#		model2d_plots(Nt, fitsum, L, 8)
#		model2d_plots(Nt, fitsum, L, 9)
#	})

	# save fit
	cat("Saving fitsum...\n")
	knots <- quantile(d, seq(0,1,length=L))
	save(L, fitsum, ud, knots, ll_calcs, DIC, pD, LPML, GG, CRPS, file=paste0("fitsums/fitsumL",L,".RData"))

	list(L=L, fitsum=fitsum, dic=dic, pD=pD)
}


})

}
