# shows that Q = M D M' in one dimension
require(fields)
require(maps)
library(multicore)
library(rstan); set_cppo("fast")
source("R/plot_1d.R")

DCT1D<-function(k){  # intrinsic
  #Spectral rep of an AR1 inverse covariance
  #t(M)%*%M = M%*%t(M) = diag(k)
  #Precision = t(M)%*%diag(d)%*%M

  M<-diag(k)
  for(j in 1:k){
     M[1,j]<-1/sqrt(k)
     for(i in 2:k){
       M[i,j]<-sqrt(2/k)*cos(pi*(i-1)*(j-1/2)/k)
  }}
  i<-1:k
  d<-2*(1-cos(pi*(i-1)/k))
  cycles<-(i-1)/2
  period<-2*k/(i-1)
list(M=M,d=d,cycles=cycles,period=period)}

"dct" <- function(n, y) {
	a <- c(1/sqrt(n), rep(sqrt(2/n), n-1))
	s <- t <- 0:(n-1)

	z <- a * sapply(t, function(time) {
		sum( cos(time * pi * (s+0.5)/n)*y )
	})
}

if (FALSE) { # verify we can decompose Q

icar <- DCT1D(4)
print(icar$M)
print(icar$d)
print(Q<-round(t(icar$M) %*% diag(icar$d) %*% icar$M,4))

print(round((icar$M) %*% Q %*% t(icar$M),4))

Q[Q == -1] <- -0.25
print(round(Q,4))
print(round((icar$M) %*% Q %*% t(icar$M),4))

done

}

if (FALSE) { # compare MY vs DCT

# simulate AR(1) time series for 100 points
set.seed(311)
n <- 100
icar <- DCT1D(n)
y <- arima.sim(n, model=list(ar = c(0.25)), sd = 1)
z1 <- icar$M %*% y
z2 <- dct(n, y)
cat("Mean difference between M x Y and DCT:", mean(z1-z2),"\n")

}

if (TRUE) { # work with data

# let's get a column of summer temps to play with (these are west to east; rows are south to north)
index <- 50
WE    <- TRUE   # use west to east or south to north?
if (WE) {
	temps <- SumTemp[,index,]
	tlon  <- lon[,index]
	tlat  <- lat[,index]
} else {
	temps <- SumTemp[index,,]
	tlon  <- lon[index,]
	tlat  <- lat[index,]
}
n     <- length(tlon)

# plot this data
if (FALSE) {
pdf("pdf/data.pdf")
	par(bty="l")
  plot(1:n,temps[,9],xlab="Grid Cell",ylab="Temp",main="Data vs Models",
		type="l",lwd=1,col="black",ylim=c(min(temps),max(temps)))   # data
  for (i in 1:1) lines(1:n,temps[,i],lwd=2,col="blue")          # GCMs
  for (i in 2:7) lines(1:n,temps[,i],lwd=0.5,col="red")         # RCMs
  for (i in 8:8) lines(1:n,temps[,i],lwd=1,col="black")         # other data
graphics.off()}

# get Zs
icar <- DCT1D(length(tlon))
z <- icar$M %*% temps

# plot Zs
if (FALSE) {
pdf("pdf/data_trans.pdf")
	par(bty="l")
  plot(icar$d,z[,9],xlab="Grid Cell",ylab="d",main="Data",type="l",lwd=2,col="black")   # data
  lines(icar$d,z[,1],lwd=2,col="blue")   # GCM
  lines(icar$d,z[,5],lwd=1,col="red")    # RCM
graphics.off()}

# plot correlations over d
d <- sapply(0:4, function(i) { mean(icar$d[1:20+i*20]) })
corrs <- sapply(0:4, function(i) {
	rev(cor(z[i*20+1:20,],z[i*20+1:20,])[,9])
})
if (FALSE) {pdf("pdf/corr.pdf")
	par(bty="l")
	plot(d, corrs[2,], xlab="d", ylab="correlation", type="p", pch=20, cex=1.0, col="black", ylim=c(-0.5,1))   # other data
	for (i in 3:8) points(d, corrs[i,], pch=20, cex=1.0, col="red")   # RCMs
	for (i in 9:9) points(d, corrs[i,], pch=20, cex=1.0, col="blue")  # GCM
graphics.off()}

# smooth correlations over d
zstar <- sqrt(icar$d)*z;
d <- sapply(0:(length(icar$d)-20), function(i) { mean(icar$d[1:20+i]) })
corrs <- sapply(0:(length(icar$d)-20), function(i) {
	rev(cor(zstar[i+1:20,],zstar[i+1:20,])[,9])
})

if (FALSE) {pdf("pdf/corr_smooth.pdf")
	par(bty="l")
	plot(d, corrs[2,], xlab="d", ylab="correlation", type="l", lwd=1.0, col="black", ylim=c(-0.5,1))   # other data
	for (i in 3:8) lines(d, corrs[i,], pch=20, lwd=1.0, col="red")   # RCMs
	for (i in 9:9) lines(d, corrs[i,], pch=20, lwd=2.0, col="blue")  # GCM
graphics.off()}
if (FALSE) {pdf("pdf/corr_smooth.pdf")
	par(bty="l")
	bw <- 0.05
	ks <- ksmooth(d, corrs[2,], "normal", bandwidth=bw)
	plot(ks$x, ks$y, xlab="d", ylab="correlation", type="l", lwd=1.0, col="black", ylim=c(-0.5,1))   # other data
	for (i in 3:8) { # RCMs
		ks <- ksmooth(d, corrs[i,], "normal", bandwidth=bw)
		lines(ks$x, ks$y, lwd=1.0, col="red")
	}
	for (i in 9:9) { # GCM
		ks <- ksmooth(d, corrs[i,], "normal", bandwidth=bw)
		lines(ks$x, ks$y, lwd=2.0, col="blue")
	}
graphics.off()}

# let's fit a model
dim1_code <- "
	data {
		int<lower=0>  n;         // sample size
		int<lower=0>  k;         // number of data sources
		matrix[n,k]   Z;         // data Z
		vector[n]     d;         // diag(D) from decomposition of Q

		int<lower=0>  L;         // number of cov matrices in mixture
		matrix[L,n]   weights;   // resolution weights for each T
	}
	transformed data {
		matrix[n,k] Zstar;
		vector[k]   zero;

		// normalize Z
		for (i in 1:n)
			Zstar[i] <- sqrt(d[i])*Z[i];

		// mean for the delta's
		for (i in 1:k)
			zero[i] <- 0.0;
	}
	parameters {
		//matrix[n,k]                  delta;
		//real<lower=0,upper=100>      sigma_e;
		vector<lower=0>[k]           omega[L];
		corr_matrix[k]               corrOmega[L];
	}
	model {
		matrix[k,k] Omega[L];
		matrix[k,k] Sigma_d[n];

		// construct Omega
		for (l in 1:L)
			for (i in 1:k)
				for (j in 1:k)
					Omega[l,i,j] <- omega[l,i]*omega[l,j]*corrOmega[l,i,j];

		// construct Sigma(T) = sum_l w_l(T) Omega_l
		for (t in 1:n) {
			Sigma_d[t] <- weights[1,t] * Omega[1];

			for (l in 2:L)
				Sigma_d[t] <- Sigma_d[t] + weights[l,t] * Omega[l];
		}

		// priors
		//sigma_e ~ uniform(0, 100);

		//for (l in 1:L)
		//	omega[l] ~ uniform(0, 100);

		for (l in 1:L)
			omega[l] ~ cauchy(0, 5);

		//for (l in 1:L)
		//	corrOmega[l] ~ lkj_corr(2.0);

		//for (t in 2:n)
		//	delta[t]' ~ multi_normal(zero, Sigma_d[t]/d[t]);

		// model
		//for (t in 1:n)
		//	Z[t] ~ normal(delta[t], sigma_e);

		for (t in 2:n)
			Zstar[t]' ~ multi_normal(zero, Sigma_d[t]);
	}

	generated quantities {
		matrix[k,k] Omega[L];
		matrix[k,k] Sigma_d[n];
		matrix[k,k] corrSigma_d[n];
		int c;

		// construct Omega
		for (l in 1:L)
			for (i in 1:k)
				for (j in 1:k)
					Omega[l,i,j] <- omega[l,i]*omega[l,j]*corrOmega[l,i,j];

		// construct Sigma(T) = sum_l w_l(T) Omega_l
		for (t in 1:n) {
			Sigma_d[t] <- weights[1,t] * Omega[1];

			for (l in 2:L)
				Sigma_d[t] <- Sigma_d[t] + weights[l,t] * Omega[l];
		}

		// construct corrSigma(T)
		for (t in 1:n)
			for (i in 1:k)
				for (j in 1:k)
						corrSigma_d[t,i,j] <- Sigma_d[t,i,j]/sqrt(Sigma_d[t,i,i]*Sigma_d[t,j,j]);
	}
"

k <- ncol(z)

# compile once...
dat = list(n=n, k=k, Z=z, d=icar$d, L=1, weights=matrix(1, nrow=1, ncol=n))
#fit1 <- stan(model_code = dim1_code, data = dat, iter = 10, chains = 1);save(fit1, file="fit1.RData");done
load("fit1.RData")

"get_weights" <- function(n, L) {
	# construct weights
	if (L > 1) {
		knots <- seq(1,n,length=L)
		sigma <- n/L

		diffs <- sapply(1:n, function(i) {
			i-knots
		})

		weights <- t(apply(diffs, 1, function(row) {
			exp(-row^2/(2*sigma^2))
		}))

		weights <- apply(weights, 2, function(col) { col/sum(col) })
	} else {
		weights <- matrix(1, nrow=L, ncol=n)
	}

	weights
}

#Ls <- 3 #1:6
Ls <- THE_L
fits <- lapply(Ls, function(L) {
	dat = list(n=n, k=k, Z=z, d=icar$d,
		L=L, weights=get_weights(n, L)
	)

	# run in parallel
	sflist <- mclapply(1:3, mc.cores=3,
		function(i) {
			stan(fit=fit1, data=dat, iter=1000,
			     chains = 1, seed=311, chain_id=i, refresh=10, verbose=FALSE,
			     #pars=c("sigma_e","omega","corrOmega","corrSigma_d"))
			     pars=c("omega","corrSigma_d")
			)
		})

	#fit <- sflist2stanfit(sflist)
	fitsum <- summary(sflist2stanfit(sflist))

	# make plots
	model1d_plots(fitsum, L, 8)
	model1d_plots(fitsum, L, 9)

	list(L=L, fitsum=fitsum)
})

}
