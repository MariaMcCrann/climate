if (!("cdata" %in% ls())) {
	source("R/load_data2.R")
}

# smooth correlations
"smooth_corr" <- function() {
	inc <- 0.005
	knots <- seq(min(f), max(f), len=(max(f)-min(f))/inc)

	cor <- vector("list", 9)
	for (i in 1:9) {
		cor[[i]] <- sapply(1:length(knots), function(k) { cor(z[f >= (knots[k]-inc)&f <= (knots[k]+inc),])[i,1:9] })
	}

	list(knots=knots, cor=cor)
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
"ls_estimates" <- function(sc, L) {
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

# smooth correlations
sc <- smooth_corr()

if (TRUE) {
	# test least squares estimates
	ests <- ls_estimates(sc, 5)
}
