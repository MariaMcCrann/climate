# simulate fake data using covariance over time
library(fda)

source("R/hmc.R")
source("R/spline_cov.R")

set.seed(03172000)
n <- 7000
t <- seq(0, 1, len=n)
k <- 9
n.off <- k*(k-1)/2

# number of basis functions
L <- 20

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
prior.theta_sd <- 100

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

grad_U <- function(x, row) { # vector of partials of U wrt x
	theta.d <- matrix(x[1:(L*k)], nrow=k, ncol=L)
	theta.o <- matrix(x[L*k+1:(L*n.off)], nrow=n.off, ncol=L)

	d_theta.d <- matrix(0, nrow=k, ncol=L)
	d_theta.o <- matrix(0, nrow=n.off, ncol=L)

	if (missing(row)) {
		# normal gradient
		rows <- 1:n
		d_theta.d <- d_theta.d -theta.d/prior.theta_sd^2
		d_theta.o <- d_theta.o -theta.o/prior.theta_sd^2
	} else {
		# only do one row and skip prior piece
		rows <- row
	}

	#for (i in 1:n) {
	for (i in rows) {
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

				P <- t(R_i) %*% D + D %*% R_i
				#T <- D %*% R_i
				#P <- t(T) + T
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

if (FALSE) { #i==12) {
print(round(t(D) %*% R_i,2))
print(round(t(R_i) %*% D,2))
}

					P <- t(R_i) %*% D + t(D) %*% R_i
					#T <- D %*% R_i
					#P <- t(T) + T
#if (i==11&&m1==1&&m2==3) { print(T) }

					#d_theta.o[index,p] <- d_theta.o[index,p] + weights[i,p]*invL_i_x_y_i[m1]*invL_i_x_y_i[m2]
#					d_theta.o[m,p] <- d_theta.o[m,p] + (t(y[i,]) %*% invSigma %*% (
#						L_i %*% t(C) + C %*% t(L_i)
#						) %*% invSigma %*% y[i,])
##					d_theta.o[index,p] <- d_theta.o[index,p] + 0.5 * t(y[i,]) %*% invSigma %*% P %*% invSigma %*% y[i,]
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

M <- function(x) { # cholesky of mass matrix M
	#chol(chol2inv(
		chol(Reduce('+', lapply(1:n, function(row) { tcrossprod(grad_U(x,row)) })))
	#))
}

s <- function(x) { # variable scales
	A <- Reduce('+', lapply(1:n, function(row) { tcrossprod(grad_U(x,row)) })) + diag(1/prior.theta_sd^2, L*(k+n.off))
	B <- chol2inv(chol(A))

	#return( (diag( B )) )
	return( sqrt(diag( B )) )

	theta.d <- matrix(x[1:(L*k)], nrow=k, ncol=L)
	theta.o <- matrix(x[L*k+1:(L*n.off)], nrow=n.off, ncol=L)

	s_theta.d <- matrix(0, nrow=k, ncol=L)
	s_theta.o <- matrix(0, nrow=n.off, ncol=L)

	for (i in 1:n) {
		R_i <- diag(exp( sapply(1:k, function(row) { crossprod(weights[i,], theta.d[row,]) }) ))
		R_i[upper.tri(R_i)] <- sapply(1:n.off, function(row) { crossprod(weights[i,], theta.o[row,]) })
#if (i==12) { print(round(R_i,2)) }

		invR_i   <- backsolve(R_i, diag(k))
		invSigma <- invR_i %*% t(invR_i)

		for (p in 1:L) {
			if (weights[i,p] == 0) next;

			for (m in 1:k) {
				D <- matrix(0, nrow=k, ncol=k)
				#D[m,m] <- exp(crossprod(weights[i,], theta.d[m,]))*weights[i,p]
				D[m,m] <- R_i[m,m]*weights[i,p]

				P <- t(R_i) %*% D + D %*% R_i

				#s_theta.d[m,p] <- s_theta.d[m,p] + (exp(crossprod(weights[i,], theta.d[m,]))*weights[i,p])^2
				#s_theta.d[m,p] <- s_theta.d[m,p] + (R_i[m,m]*weights[i,p])^2
				s_theta.d[m,p] <- s_theta.d[m,p] + sum(diag( invSigma %*% P %*% invSigma %*% P ))
			}

			index <- 1
			for (m1 in 1:(k-1)) {
				for (m2 in (m1+1):k) {
					D <- matrix(0, nrow=k, ncol=k)
					D[m1,m2] <- weights[i,p]

					P <- t(R_i) %*% D + t(D) %*% R_i
#if (i==2&&index==1) { print(round(invSigma %*% P, 2)) }
#if (index==1) { print(round(diag( invSigma %*% P %*% invSigma %*% P ),2)) }
#if (i==2&&index==1) { cat(round( (diag( invSigma %*% P %*% invSigma %*% P )), 2),"\n") }
					s_theta.o[index,p] <- s_theta.o[index,p] + sum(diag( invSigma %*% P %*% invSigma %*% P ))

					index <- index+1
				}
			}
 

		}
#cat("\n")
#print(c(i,round(s_theta.o[1],2)))
	}

#print(round(as.vector(s_theta.o),2))

	s_theta.d <- 0.5*s_theta.d + 1/prior.theta_sd^2
	s_theta.o <- 0.5*s_theta.o + 1/prior.theta_sd^2

	sds <- 1/sqrt((c(as.vector(s_theta.d), as.vector(s_theta.o))))
	#sds <- max(sds)/sds

	#sqrt(1/(c(as.vector(s_theta.d), as.vector(s_theta.o))))
	#1/sds

	#sds/2
	#sqrt((c(as.vector(s_theta.d), as.vector(s_theta.o))))
}

v1 <- c(as.vector(alpha.d), as.vector(alpha.o))
v2 <- rnorm(L*(k+n.off))
v3 <- rep(0, L*(k+n.off))
#print( c(U(v1), U(v2), U(v3)) )
#print( round(grad_U(v1),3) ); done; #print( round(grad_U(v2),3) ); print( round(grad_U(v3),3) );done
#print( round(s(v1),3) ); print( round(s(v2),3) ); print( round(s(v3),3) );done
#( round(s(v1),3) ); ( round(s(v2),3) ); ( round(s(v3),3) );done

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
sapply(1:(L*(k+n.off)), function(w) {
	cat("w=",w,"\n")
	e <- rep(0, length(v1))
	e[w] <- 0.0001
	print( c( (U(v1+e)-U(v1))/e[w], grad_U(v1)[w]) )
	print( c( (U(v2+e)-U(v2))/e[w], grad_U(v2)[w]) )
	print( c( (U(v3+e)-U(v3))/e[w], grad_U(v3)[w]) )
})
done
}

if (TRUE) {
# construct weights
if (L == 1) {
	weights <- matrix(1, nrow=n, ncol=1)
} else {
	basis   <- create.bspline.basis(c(min(t),max(t)),norder=4,nbasis=L)
	knots   <- knots(basis)
	weights <- getbasismatrix(t, basis, nderiv=0)
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

t1 <- proc.time()
for (i in 1:1) print(round(U(v3),2))
#for (i in 1:1) print(round(grad_U(v1),2))
#for (i in 1:1) print(round(s(v1),2))
#for (i in 1:50) { U(v1); grad_U(v1); s(v1); }
print(proc.time()-t1)

t1 <- proc.time()
print(spline_cov_lk(data=list(prior=prior.theta_sd,
		n=n, k=k, y=y, L=L, Nnz=Nnz, Mnz=Mnz-1, Wnz=Wnz), eval=v3)$lk)
print(proc.time()-t1)
done
}

if (TRUE) {

"sim_fit" <- function(fitL, step_e, step_L) {
	# construct weights
	if (fitL == 1) {
		weights <- matrix(1, nrow=n, ncol=1)
	} else {
		basis   <- create.bspline.basis(c(min(t),max(t)),norder=4,nbasis=fitL)
		knots   <- knots(basis)
		weights <- getbasismatrix(t, basis, nderiv=0)
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

	Niter <- 100000

	samples <- matrix(0, nrow=Niter, ncol=fitL*(k+n.off))

	inits <- v1; #rep(0, fitL*(k+n.off))
	t1 <- proc.time()
	fit <- spline_cov(data=list(prior=prior.theta_sd,
		n=n, k=k, y=y, L=fitL, Nnz=Nnz, Mnz=Mnz-1, Wnz=Wnz),
	  step_e=step_e, step_L=step_L, inits=inits, Niter=Niter, thin=1, verbose=TRUE)
	print(proc.time()-t1)

	#dic <- spline_cov_dic(fit, 100)
	res <- mcmc( matrix(fit$samples, nrow=Niter) )

	list(fit=fit, res=res) #, dic=dic)
}

set.seed(1983)
fit <- sim_fit(L, 0.20, 10)
#fit1 <- sim_fit(5, 0.1, 25)
#fit2 <- sim_fit(10, 0.1, 25)
#fit3 <- sim_fit(15, 0.1, 25)
#fit4 <- sim_fit(20, 0.1, 25)

#} else {

start <- v1
step.e <- 0.05
step.L <- 10 #round(1/step.e)
print(c(step.e, step.L, step.L*step.e))

set.seed(1983)
nsamples <- 5
samples <- matrix(0, nrow=nsamples, ncol=L*(k+n.off))
current <- start
acc <- 0
for (i in 1:nsamples) {
	res <- HMC(U, grad_U, step.e, step.L, current, s=s)
	#res <- HMC(U, grad_U, step.e, step.L, current, M=M)
	current <- res$val
#print(c(res$acc,current))
	samples[i,] <- current
	acc <- acc+res$acc
}
done
mcmc.samples <- mcmc(samples) #samples[round(nsamples/2):nsamples,])
cat("HMC acceptance rate:");print(acc/nsamples)
print(summary(mcmc.samples))
#print(acf(mcmc.samples,plot=FALSE))
cat("Effective sample size:\n");print(effectiveSize(mcmc.samples))

}
