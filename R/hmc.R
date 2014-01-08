# basic HMC code adapted from page 125 of "Handbook of Markov Chain Monte Carlo"
require(coda)
require(MASS)

"HMC" <- function(U, grad_U, epsilon, L, current_q, s) {
	q <- current_q
	p <- rnorm(length(q), 0, 1)  # independent std normal variates
	current_p <- p

	if (!missing(s)) { # different step sizes
		epsilon <- epsilon*s(q)
#print(epsilon);#done
	}

	# make a half step for momentum at the beginning
	p <- p - epsilon * grad_U(q) / 2

	# alternate full steps for position and momentum
	for (i in 1:L) {
		# make a full step for the position
		q <- q + epsilon * p
		# make a full step for the momentum, except at the end of trajectory
		if (i != L) {
			p <- p - epsilon * grad_U(q)
		}
	}

	# make a half step for momentum at the end
	p <- p - epsilon * grad_U(q) / 2
	# negate momentum at end of trajectory to make the proposal symmetric
	p <- -p

#print(tail(current_q)); print(tail(q))
#print(current_q); print(q);

	# evaluate potentail and kinetic energies at start and end of trajectory
	proposed_U <- U(q)
	proposed_K <- sum(p^2)/2
	current_U <- U(current_q)
	current_K <- sum(current_p^2)/2

	# accept or reject the state at end of trajectory, returning either
	# the position at the end of the trajectory or the initial position
	pr <- log(runif(1))
	vs <- current_U-proposed_U + current_K-proposed_K

if (TRUE) {
cat("Current:");print(c(current_U,current_K))
cat("Proposed:");print(c(proposed_U,proposed_K))
print(round(epsilon[1:10],3));
print(round(current_q[1:10],3)); print(round(q[1:10],3));
#print(round(p^2,3));
#print(c(pr=round(pr,3),vs=round(vs,3)))
print(c(pr,vs,pr<vs));
#done
}

if (is.nan(vs)) {
	print(c(current_U,proposed_U,current_K,proposed_K))
#	print(U)
#	print(grad_U)
#	print(c(epsilon,L,current_q))
}

	if (!is.nan(vs) && pr < vs) {
		return(list(value=q,p=p,acc=1)) # accept
	} else {
		return(list(value=current_q,p=current_p,acc=0))  # reject
	}
}

if (0) { # simple multinomial

# generate and fit sample data
set.seed(311)
mu <- c(0,-1.0,-0.1,-0.5)
probs <- exp(mu)/sum(exp(mu))

n <- 5000
y <- as.vector(rmultinom(1,n,probs))

print(probs)
print(y)

"U" <- function(x) { # -log kernel(x)
	# kernel is [ 1 + sum_l exp(mu_l) ]^(-n) prod_i>1 exp(mu_i)^y_i
	# log k = -n log[1 + sum_l exp(mu_l)] + sum_i>1 y_i * mu_i
	# -log k = n log[ 1 + sum_l exp(mu_l) ] - sum_i>1 y_i*mu_i
	n * log(1 + sum(exp(x))) - sum(y[-1] * x)
}

"grad_U" <- function(x) { # vector of partials of U
	# partial wrt mu_i: n * exp(mu_i) / [ 1 + sum_l exp(mu_l) ] - y_i
	expx <- exp(x)
	denom <- 1 + sum(expx)
	n * expx / denom - y[-1]
}

# fit with HMC
#epsilon <- 0.005; L <- 25  ### gives negative correlation!!!
epsilon <- 0.005; L <- 25  ### gives negative correlation!!!
print(L*epsilon)

nsamples <- 25000
samples <- matrix(0, nrow=nsamples, ncol=3)
current <- c(0,0,0)
acc <- 0
for (i in 1:nsamples) {
	res <- HMC(U, grad_U, epsilon, L, current)
	current <- res$val
	samples[i,] <- current
	acc <- acc+res$acc
}
mcmc.samples <- mcmc(samples)
cat("Acceptance rate:");print(acc/nsamples)
print(summary(mcmc.samples))
acf(mcmc.samples)
print(effectiveSize(mcmc.samples))

}

if (0) {  # the bivariate normal example from the book

# sample from bivariate normal with mean 0, std 1, and corr=0.95
invSigma <- chol2inv(chol( matrix( c(1,0.95,0.95,1) ,nrow=2) ))
U <- function(x) {  # -log K(x)
	0.5 * ( t(x) %*% invSigma %*% x )
}

grad_U <- function(x) { # vector of partials of U wrt X
	# U looks like:
	# -0.5 * [ x_1*s_11 + x_2*s_12, x_1*s_12 + x_2*s_22 ] x
	# = -0.5 * [ x_1^2*s_11 + x_1*x_2*s_12 + x_1*x_2*s_12 + x_2^2*s_22 ]
	# = -0.5 * [ x_1^2*s_11 + 2*x_1*x_2*s_12 + x_2^2*s_22 ]
	c(x[1]*invSigma[1,1] +x[2]*invSigma[1,2], x[2]*invSigma[2,2] +x[1]*invSigma[1,2])
}


start <- c(0,0)
epsilon <- 0.25
L <- 25
print(L*epsilon)

set.seed(311)
nsamples <- 25000
samples <- matrix(0, nrow=nsamples, ncol=2)
current <- start
acc <- 0
for (i in 1:nsamples) {
	res <- HMC(U, grad_U, epsilon, L, current)
	current <- res$val
	samples[i,] <- current
	acc <- acc+res$acc
}
mcmc.samples <- mcmc(samples)
cat("HMC acceptance rate:");print(acc/nsamples)
print(summary(mcmc.samples))
acf(mcmc.samples)
cat("Effective sample size:\n");print(effectiveSize(mcmc.samples))

# compare to MH
acc <- 0
"MH" <- function(k, loc) {
	new <- mvrnorm(1, mu=loc, Sigma=0.5*diag(2));
	pr <- log(runif(1))
	vs <- -U(new) - -U(loc)
	if (pr < vs) {
		acc <<- acc+1
		return(new);
	} else {
		return(loc);
	}
}

mh.samples <- matrix(0, nrow=nsamples, ncol=2)
current <- start
for (i in 1:nsamples) {
	current <- MH(U, current)
	mh.samples[i,] <- current
}
mcmc.mh.samples <- mcmc(mh.samples)
cat("M-H acceptance rate:");print(acc/nsamples)
print(summary(mcmc.mh.samples))
cat("Effective sample size:\n");print(effectiveSize(mcmc.mh.samples))

}

if (0) {  # unknown mean/var from pg 64 of BDA 3

# sample data
set.seed(1983)
data   <- list()
data$n    <- 100
data$y    <- rnorm(data$n, mean=10, sd=3)
data$ybar <- mean(data$y)
data$s2   <- var(data$y)

U <- function(x) {  # -log kernel(x)
	#lk <- -(data$n+1)*x[2] -0.5*exp(-2*x[2])*( (data$n-1)*data$s2 + data$n*(data$ybar - x[1])^2 )
	lk <- -data$n*x[2] -0.5*exp(-2*x[2])*( (data$n-1)*data$s2 + data$n*(data$ybar - x[1])^2 )

	-lk
}

grad_U <- function(x) { # vector of partials of U wrt x
	-c(
		data$n*exp(-2*x[2])*(data$ybar-x[1]),
		#-(data$n+1)+exp(-2*x[2])*( (data$n-1)*data$s2 + data$n*(data$ybar-x[1])^2 )
		-data$n+exp(-2*x[2])*( (data$n-1)*data$s2 + data$n*(data$ybar-x[1])^2 )
	)
}

s <- function(x) { # variable scales
	sqrt( 1/c(
		data$n*exp(-2*x[2]),
		#2*data$n*(1+x[1]^2*exp(-2*x[2]))
		#2*data$n*exp(-2*x[2])
		2*data$n
	) )
}

start <- c(mean(data$y),log(sd(data$y)))
epsilon <- .25
L <- round(2/epsilon)
print(L*epsilon)

set.seed(311)
nsamples <- 1000
samples <- matrix(0, nrow=nsamples, ncol=2)
current <- start
acc <- 0
for (i in 1:nsamples) {
	res <- HMC(U, grad_U, epsilon, L, current, s=s)
	current <- res$val
	samples[i,] <- current
	acc <- acc+res$acc
}
mcmc.samples <- mcmc(samples)
cat("HMC acceptance rate:");print(acc/nsamples)
print(summary(mcmc.samples))
acf(mcmc.samples)
cat("Effective sample size:\n");print(effectiveSize(mcmc.samples))

done

# compare to MH
acc <- 0
"MH" <- function(k, loc) {
	new <- mvrnorm(1, mu=loc, Sigma=0.5*diag(2));
	pr <- log(runif(1))
	vs <- -U(new) - -U(loc)
	if (pr < vs) {
		acc <<- acc+1
		return(new);
	} else {
		return(loc);
	}
}

mh.samples <- matrix(0, nrow=nsamples, ncol=2)
current <- start
for (i in 1:nsamples) {
	current <- MH(U, current)
	mh.samples[i,] <- current
}
mcmc.mh.samples <- mcmc(mh.samples)
cat("M-H acceptance rate:");print(acc/nsamples)
print(summary(mcmc.mh.samples))
cat("Effective sample size:\n");print(effectiveSize(mcmc.mh.samples))

}
