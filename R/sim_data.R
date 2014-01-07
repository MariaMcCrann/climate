# simulate fake data using covariance over time
library(fda)

source("R/hmc.R")

set.seed(03172000)
n <- 100
t <- seq(0, 1, len=n)
k <- 3
n.off <- k*(k-1)/2

# number of basis functions
L <- 5

if (L > 1) {
	# construct weights
	basis   <- create.bspline.basis(c(min(t),max(t)),norder=4,nbasis=L)
	knots   <- knots(basis)
	weights <- getbasismatrix(t, basis, nderiv=0)
} else {
	weights <- matrix(1, nrow=n, ncol=1)
}

# generate coefficients
alpha.d <- matrix(rnorm(L*k), nrow=k, ncol=L)
alpha.o <- matrix(rnorm(L*n.off), nrow=n.off, ncol=L)

# simulate data
y <- t( sapply(1:n, function(i) {
#	L_i <- diag(exp( sapply(1:k, function(row) { crossprod(weights[i,], alpha.d[row,]) }) ))
#	L_i[lower.tri(L_i)] <- sapply(1:n.off, function(row) { crossprod(weights[i,], alpha.o[row,]) })
#	Sigma_i <- L_i %*% t(L_i)
#	L_i %*% rnorm(k)

	R_i <- diag(exp( sapply(1:k, function(row) { crossprod(weights[i,], alpha.d[row,]) }) ))
	R_i[upper.tri(R_i)] <- sapply(1:n.off, function(row) { crossprod(weights[i,], alpha.o[row,]) })
	t(R_i) %*% rnorm(k)
}) )

# prior SD
prior.theta_sd <- 1

U.old <- function(x) {  # -log kernel(x)
	theta.d <- matrix(x[1:(L*k)], nrow=k, ncol=L)
	theta.o <- matrix(x[L*k+1:(L*n.off)], nrow=n.off, ncol=L)

	lk <- -0.5*sum(x^2)/prior.theta_sd^2 + # prior
		sum(sapply(1:n, function(i) {
			L_i <- diag(exp( sapply(1:k, function(row) { crossprod(weights[i,], theta.d[row,]) }) ))
			L_i[lower.tri(L_i)] <- sapply(1:n.off, function(row) { crossprod(weights[i,], theta.o[row,]) })

			-sapply(1:k, function(row){ crossprod(weights[i,], theta.d[row,]) }) -0.5*t(y[i,]) %*% chol2inv(t(L_i)) %*% y[i,]
		 }))

	-lk
}

U <- function(x) {  # -log kernel(x)
	theta.d <- matrix(x[1:(L*k)], nrow=k, ncol=L)
	theta.o <- matrix(x[L*k+1:(L*n.off)], nrow=n.off, ncol=L)

	cstd <- colSums(theta.d)

	lk <- -0.5*sum(x^2)/prior.theta_sd^2 + # prior
		sum(sapply(1:n, function(i) {
			R_i <- diag(exp( sapply(1:k, function(row) { crossprod(weights[i,], theta.d[row,]) }) ))
			R_i[upper.tri(R_i)] <- sapply(1:n.off, function(row) { crossprod(weights[i,], theta.o[row,]) })

			a <- backsolve(R_i, y[i,], transpose=TRUE)

			#c <- -crossprod(weights[i,], cstd) -0.5*t(y[i,]) %*% chol2inv(R_i) %*% y[i,]
			c <- -crossprod(weights[i,], cstd) -0.5* crossprod(a)
		 }))

	-lk
}

grad_U.old <- function(x) { # vector of partials of U wrt x
	theta.d <- matrix(x[1:(L*k)], nrow=k, ncol=L)
	theta.o <- matrix(x[L*k+1:(L*n.off)], nrow=n.off, ncol=L)

	d_theta.d <- matrix(0, nrow=k, ncol=L)
	d_theta.o <- matrix(0, nrow=n.off, ncol=L)

	for (i in 1:n) {
		L_i <- diag(exp( sapply(1:k, function(row) { crossprod(weights[i,], theta.d[row,]) }) ))
		L_i[lower.tri(L_i)] <- sapply(1:n.off, function(row) { crossprod(weights[i,], theta.o[row,]) })

		#invL_i_x_y_i <- forwardsolve(L_i, y[i,])
		invL_i       <- forwardsolve(L_i, diag(k))
		invL_i_x_y_i <- invL_i %*% y[i,]

		invSigma <- t(invL_i) %*% invL_i

		F <- t(invL_i) %*% invL_i_x_y_i

		for (p in 1:L) {
			for (m in 1:k) {
				d_theta.d[m,p] <- d_theta.d[m,p] -weights[i,p]

				if (weights[i,p] > 0) {
					a_ip <- exp(crossprod(weights[i,], theta.d[m,]))*weights[i,p]

					C <- matrix(0, nrow=k, ncol=k)
					C[m,m] <- a_ip
					#d_theta.d[m,p] <- d_theta.d[m,p] -weights[i,p] + t(invL_i_x_y_i) %*% C %*% t(invL_i) %*% invL_i_x_y_i
					d_theta.d[m,p] <- d_theta.d[m,p] + (t(y[i,]) %*% invSigma %*% (
						L_i %*% t(C) + C %*% t(L_i)
						) %*% invSigma %*% y[i,])
#print(L_i %*% t(C))
#print(C %*% t(L_i))
				}
			}

			index <- 1
			for (m1 in 2:k) {
				for (m2 in 1:(m1-1)) {
					C <- matrix(0, nrow=k, ncol=k)
					C[m1,m2] <- C[m2,m1] <- weights[i,p]
					#d_theta.o[index,p] <- d_theta.o[index,p] + weights[i,p]*invL_i_x_y_i[m1]*invL_i_x_y_i[m2]
					d_theta.o[m,p] <- d_theta.o[m,p] + (t(y[i,]) %*% invSigma %*% (
						L_i %*% t(C) + C %*% t(L_i)
						) %*% invSigma %*% y[i,])

					index <- index+1
				}
			}
		}
	}

	d_theta.d <- d_theta.d -theta.d/prior.theta_sd^2
	d_theta.o <- d_theta.o -theta.o/prior.theta_sd^2

   -c(as.vector(d_theta.d), as.vector(d_theta.o))
}

grad_U <- function(x) { # vector of partials of U wrt x
	theta.d <- matrix(x[1:(L*k)], nrow=k, ncol=L)
	theta.o <- matrix(x[L*k+1:(L*n.off)], nrow=n.off, ncol=L)

	d_theta.d <- matrix(0, nrow=k, ncol=L)
	d_theta.o <- matrix(0, nrow=n.off, ncol=L)

	d_theta.d <- d_theta.d -theta.d/prior.theta_sd^2
	d_theta.o <- d_theta.o -theta.o/prior.theta_sd^2

	for (i in 1:n) {
		R_i <- diag(exp( sapply(1:k, function(row) { crossprod(weights[i,], theta.d[row,]) }) ))
		R_i[upper.tri(R_i)] <- sapply(1:n.off, function(row) { crossprod(weights[i,], theta.o[row,]) })

		#invL_i_x_y_i <- forwardsolve(L_i, y[i,])
		invR_i       <- backsolve(R_i, diag(k))
		#invR_i_x_y_i <- invR_i %*% y[i,]
		A <- invR_i %*% t(invR_i)
		A_y_i <- A %*% y[i,]

if (FALSE) { #i == 12) {
print(round(R_i,2))
print(round(invR_i,2))
print(round(A,2))
print(round(A_y_i,2))
}

		invSigma <- invR_i %*% t(invR_i)

		for (p in 1:L) {
			if (weights[i,p] == 0) { next; }

			for (m in 1:k) {
				d_theta.d[m,p] <- d_theta.d[m,p] -weights[i,p]

				D <- matrix(0, nrow=k, ncol=k)
				#D[m,m] <- exp(crossprod(weights[i,], theta.d[m,]))*weights[i,p]
				D[m,m] <- R_i[m,m]*weights[i,p]

				#P <- t(R_i) %*% D + D %*% R_i
				T <- D %*% R_i
				P <- t(T) + T
if (FALSE) { #i == 11 && m==3) {
	print(c(i,p,m));
#	print(round(P,2));
	print( round(t(A_y_i) %*% P %*% A_y_i,2) );
}

					#d_theta.d[m,p] <- d_theta.d[m,p] + t(invR_i_x_y_i) %*% D %*% invR_i %*% invR_i_x_y_i
					#d_theta.d[m,p] <- d_theta.d[m,p] + 0.5 * t(invR_i_x_y_i) %*% D %*% invR_i %*% invR_i_x_y_i
##					d_theta.d[m,p] <- d_theta.d[m,p] + 0.5 * t(y[i,]) %*% invSigma %*% P %*% invSigma %*% y[i,]
				d_theta.d[m,p] <- d_theta.d[m,p] + 0.5 * t(A_y_i) %*% P %*% A_y_i
			}

			index <- 1
			for (m1 in 1:(k-1)) {
				for (m2 in (m1+1):k) {
					D <- matrix(0, nrow=k, ncol=k)
					D[m1,m2] <- weights[i,p]

					#P <- t(R_i) %*% D + t(D) %*% R_i
					T <- D %*% R_i
					P <- t(T) + T
#if (i==11&&m1==1&&m2==3) { print(T) }

					#d_theta.o[index,p] <- d_theta.o[index,p] + weights[i,p]*invL_i_x_y_i[m1]*invL_i_x_y_i[m2]
#					d_theta.o[m,p] <- d_theta.o[m,p] + (t(y[i,]) %*% invSigma %*% (
#						L_i %*% t(C) + C %*% t(L_i)
#						) %*% invSigma %*% y[i,])
					#d_theta.o[index,p] <- d_theta.o[index,p] + 0.5 * t(y[i,]) %*% invSigma %*% P %*% invSigma %*% y[i,]
					d_theta.o[index,p] <- d_theta.o[index,p] + 0.5 * t(A_y_i) %*% P %*% A_y_i

					index <- index+1
				}
			}
		}
	}

   -c(as.vector(d_theta.d), as.vector(d_theta.o))
}

s.old <- function(x) { # variable scales
	theta.d <- matrix(x[1:(L*k)], nrow=k, ncol=L)
	theta.o <- matrix(x[L*k+1:(L*n.off)], nrow=n.off, ncol=L)

	s_theta.d <- matrix(0, nrow=k, ncol=L)
	s_theta.o <- matrix(0, nrow=n.off, ncol=L)

	for (i in 1:n) {
		L_i <- diag(exp( sapply(1:k, function(row) { crossprod(weights[i,], theta.d[row,]) }) ))
		L_i[lower.tri(L_i)] <- sapply(1:n.off, function(row) { crossprod(weights[i,], theta.o[row,]) })
		#invL_i <- forwardsolve(L_i, diag(k))

		for (p in 1:L) {
			for (m in 1:k) {
				a_ip <- exp(crossprod(weights[i,], theta.d[m,]))*weights[i,p]

				s_theta.d[m,p] <- s_theta.d[m,p] +a_ip^2
			}

			index <- 1
			for (m1 in 2:k) {
				for (m2 in 1:(m1-1)) {
					s_theta.o[index,p] <- s_theta.o[index,p] + weights[i,p]^2

					index <- index+1
				}
			}
		}
	}

	s_theta.d <- 2*s_theta.d + 1/prior.theta_sd^2
	s_theta.o <-   s_theta.o + 1/prior.theta_sd^2

	sqrt(1/c(as.vector(s_theta.d), as.vector(s_theta.o)))
}

s <- function(x) { # variable scales
	theta.d <- matrix(x[1:(L*k)], nrow=k, ncol=L)
	theta.o <- matrix(x[L*k+1:(L*n.off)], nrow=n.off, ncol=L)

	s_theta.d <- matrix(0, nrow=k, ncol=L)
	s_theta.o <- matrix(0, nrow=n.off, ncol=L)

	for (i in 1:n) {
		R_i <- diag(exp( sapply(1:k, function(row) { crossprod(weights[i,], theta.d[row,]) }) ))
		R_i[upper.tri(L_i)] <- sapply(1:n.off, function(row) { crossprod(weights[i,], theta.o[row,]) })

		for (p in 1:L) {
			for (m in 1:k) {
				s_theta.d[m,p] <- s_theta.d[m,p] + (exp(crossprod(weights[i,], theta.d[m,]))*weights[i,p])^2
			}

			index <- 1
			for (m1 in 2:k) {
				for (m2 in 1:(m1-1)) {
					s_theta.o[index,p] <- s_theta.o[index,p] + weights[i,p]^2

					index <- index+1
				}
			}
		}
	}

	s_theta.d <- 2*s_theta.d + 1/prior.theta_sd^2
	s_theta.o <-   s_theta.o + 1/prior.theta_sd^2

	sds <- 1/sqrt((c(as.vector(s_theta.d), as.vector(s_theta.o))))
	#sds <- max(sds)/sds

	#sqrt(1/(c(as.vector(s_theta.d), as.vector(s_theta.o))))
	#1/sds

	sds/2
	#sqrt((c(as.vector(s_theta.d), as.vector(s_theta.o))))
}

v1 <- c(as.vector(alpha.d), as.vector(alpha.o))
v2 <- rnorm(L*(k+n.off))
v3 <- rep(0, L*(k+n.off))
#print( c(U(v1), U(v2), U(v3)) )
#print( round(grad_U(v1),3) ); done; #print( round(grad_U(v2),3) ); print( round(grad_U(v3),3) );done
#print( round(s(v1),3) ); print( round(s(v2),3) ); print( round(s(v3),3) )

#eps/sd < 2 => eps < 2*sd
#done

#print(U(v1)); print(U(v2)); print(U(v3))

#eps<-0.0001; print( c( (U(c(v1[1]+eps,v1[-1]))-U(v1))/eps, grad_U(v1)[1]) )
#eps<-0.0001; print( c( (U(c(v2[1]+eps,v2[-1]))-U(v2))/eps, grad_U(v2)[1]) )
#eps<-0.0001; print( c( (U(c(v3[1]+eps,v3[-1]))-U(v3))/eps, grad_U(v3)[1]) )

#eps<-0.0001; print( c( (U(c(v1[1:k],v1[k+1]+eps,v1[(k+2):length(v1)]))-U(v1))/eps, grad_U(v1)[k+1]) )
#eps<-0.0001; print( c( (U(c(v2[1:k],v2[k+1]+eps,v2[(k+2):length(v2)]))-U(v2))/eps, grad_U(v2)[k+1]) )
#eps<-0.0001; print( c( (U(c(v3[1:k],v3[k+1]+eps,v3[(k+2):length(v3)]))-U(v3))/eps, grad_U(v3)[k+1]) )

if (FALSE) { # check gradient
#w <- L*k+1
sapply(1:(L*(k+n.off)-2), function(w) {
	cat("w=",w,"\n")
	eps<-0.0001; print( c( (U(c(v1[1:w],v1[w+1]+eps,v1[(w+2):length(v1)]))-U(v1))/eps, grad_U(v1)[w+1]) )
	eps<-0.0001; print( c( (U(c(v2[1:w],v2[w+1]+eps,v2[(w+2):length(v2)]))-U(v2))/eps, grad_U(v2)[w+1]) )
	eps<-0.0001; print( c( (U(c(v3[1:w],v3[w+1]+eps,v3[(w+2):length(v3)]))-U(v3))/eps, grad_U(v3)[w+1]) )
})
done
}

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

Niter <- 100
samples <- matrix(0, nrow=Niter, ncol=L*(k+n.off))

t1 <- proc.time()
#for (i in 1:1) U(v1) #print(round(U(v1),2))
for (i in 1:1) print(round(grad_U(v3),2))
print(proc.time()-t1)

source("R/spline_cov.R")
t1 <- proc.time()
for (i in 1:1) fit <- spline_cov(prior=1,
	n=n, k=k, y=y, L=L, Nnz=Nnz, Mnz=Mnz, Wnz=Wnz,
  step_e=0.5, step_L=10, inits=v3, Niter=Niter, samples=samples)
print(proc.time()-t1)

done

start <- v1
step.e <- 1
step.L <- 10 #round(1/step.e)
print(c(step.e, step.L, step.L*step.e))

set.seed(311)
nsamples <- 500
samples <- matrix(0, nrow=nsamples, ncol=L*(k+n.off))
current <- start
acc <- 0
for (i in 1:nsamples) {
	res <- HMC(U, grad_U, step.e, step.L, current, s=s)
	current <- res$val
print(c(res$acc,current))
	samples[i,] <- current
	acc <- acc+res$acc
}
mcmc.samples <- mcmc(samples[round(nsamples/2):nsamples,])
cat("HMC acceptance rate:");print(acc/nsamples)
print(summary(mcmc.samples))
print(acf(mcmc.samples,plot=FALSE))
cat("Effective sample size:\n");print(effectiveSize(mcmc.samples))
