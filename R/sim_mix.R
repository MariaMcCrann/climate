# simulate fake data using covariance mixture
library(fda)
library(MCMCpack)

#source("R/mix_cov.R")

set.seed(03172000)
n <- 1000
t <- seq(0, 1, len=n)
k <- 9

# number in mixture
L <- 20

if (L > 1) {
	# construct weights
	basis   <- create.bspline.basis(c(min(t),max(t)),norder=4,nbasis=L)
	knots   <- knots(basis)
	weights <- getbasismatrix(t, basis, nderiv=0)
} else {
	weights <- matrix(1, nrow=n, ncol=1)
}

# generate covs
nu <- 15
covs <- vector("list", L)
for (i in 1:L) {
	covs[[i]] <- riwish(nu, (nu-k-1) * diag(k))
}

# simulate data
y <- t( sapply(1:n, function(i) {
	Sigma_i <- Reduce('+', mapply("*", covs, weights[i,], SIMPLIFY=FALSE))

	t(chol(Sigma_i)) %*% rnorm(k)
}) )

"fit_mix" <- function(data, prior, Niter=500) {
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
	keep.alpha  <- array(NA, dim=c(Niter,data$n,data$L,data$k))
	keep.ll     <- rep(NA, Niter)

	# means of each observation
	mu <- matrix(0, nrow=data$n, ncol=data$k)

	# compute mean of each observation
	sapply(1:data$n, function(i) {
		mu[i,] <<- 0
		sapply(which(data$w[i,] > 0), function(l) {
			if (l != m_w[i]) {
				mu[i,] <<- mu[i,] + rt_w[i,l] * alpha[i,l,]
			}
		})
	})

	for (iter in 1:Niter) {
if (TRUE) {
		# update Omega
		sapply(1:L, function(l) {
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
			B       <- prior$S + S_y + S_alpha
if (l==1) {
	#print( a-data$k-1 )
	#print(round(S_alpha,3))
	#print(a)
	print( round(B/(a-data$k-1),3) )
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

		# compute mean of each observation
		sapply(1:data$n, function(i) {
			mu[i,] <<- 0
			sapply(which(data$w[i,] > 0), function(l) {
				if (l != m_w[i]) {
					mu[i,] <<- mu[i,] + rt_w[i,l] * alpha[i,l,]
				}
			})
		})

		# compute log-lik
if (FALSE) {
		ll <- sum(sapply(1:data$n, function(i) {
			invS <- invOmega[m_w[i],,]/data$w[i,m_w[i]]
			cholInvS <- chol(invS)
#print(log(det( solve(invS) )))
#print(log(1/det(invS)))
#print(2*sum(log(diag(chol(solve(invS))))))
#print(2*sum(-log(diag(cholInvS))))
			-(sum(-log(diag(cholInvS)))) -0.5*t(data$y[i,]-mu[i,]) %*% invS %*% (data$y[i,]-mu[i,])
		}))
} else {
		ll <- sum(sapply(1:data$n, function(i) {
			S <- Reduce('+', lapply(which(data$w[i,] > 0), function(l) { data$w[i,l]*Omega[l,,] }))
			cholS <- chol(S)
			invS <- chol2inv(cholS)
			-sum(log(diag(cholS))) -0.5*t(data$y[i,]) %*% invS %*% (data$y[i,])
		}))
}

		if (iter %% 10==0) print(iter)

		keep.Omega[iter,,,] <- Omega
		keep.alpha[iter,,,] <- alpha
		keep.ll[iter] <- ll
	}

	list(Omega=keep.Omega, alpha=keep.alpha, ll=keep.ll)
}

set.seed(311)
fit <- fit_mix(data=list(y=y, n=n, k=k, L=L, w=weights),
	prior=list(nu=15, S=(15-k-1)*diag(k))
)
