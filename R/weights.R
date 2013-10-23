# test out weight functions

"w_gauss" <- function(d, L, knots) {
	# construct weights from Gaussian kernels
	if (L > 1) {
		sigma <- 0.5*max(d)/L

		diffs <- sapply(d, function(i) {
			i-knots
		})

		weights <- t(apply(diffs, 1, function(row) {
			exp(-row^2/(2*sigma^2))
		}))

		weights <- apply(weights, 2, function(col) { col/sum(col) })
	}

	t(weights)
}

if (FALSE) {
	L <- 4
	#knots <- seq(min(d), max(d), len=L)
	knots <- c(min(d), 0.1, max(d)/2, max(d))
	ud <- seq(min(d), max(d), len=100)
	w.gauss <- w_gauss(ud, L, knots)
	plot(ud, w.gauss[,1], type="l", ylim=c(0,1));
	for(j in 2:L) { lines(ud,w.gauss[,j],lty=j) }
}

"w_bgauss" <- function(d, L, knots) {
	# construct weights from truncated Gaussian kernels
	if (L > 1) {
		sigma <- 0.5*max(d)/L

		diffs <- sapply(d, function(i) {
			i-knots
		})

		# create bounds
		bounds.l <- rep(0, length(knots))
		bounds.u <- rep(0, length(knots))

		for (i in 1:length(knots)) {
			if (i == 1) {
				bounds.l[i] <- min(d)
				bounds.u[i] <- 0.5*(knots[i]+knots[i+1])
			} else if (i == length(knots)) {
				bounds.l[i] <- 0.5*(knots[i-1]+knots[i])
				bounds.u[i] <- max(d)
			} else {
				bounds.l[i] <- 0.5*(knots[i-1]+knots[i])
				bounds.u[i] <- 0.5*(knots[i]+knots[i+1])
			}
		}

		in_lower <- t(sapply(d, function(i) { as.integer(i >= bounds.l) }))
		in_upper <- t(sapply(d, function(i) { as.integer(i <= bounds.u) }))

#print(head(in_lower)); print(head(in_upper))
#print(rbind(bounds.l,bounds.u));done

		weights <- t(sapply(1:nrow(diffs), function(row) {
			exp(-diffs[row,]^2/(2*sigma^2)) * in_lower[,row] * in_upper[,row]
		}))

		weights <- apply(weights, 2, function(col) { col/sum(col) })
	}

	t(weights)
}

if (FALSE) {
	L <- 4
	#knots <- seq(min(d), max(d), len=L)
	knots <- c(min(d), 0.1, max(d)/2, max(d))
	ud <- seq(min(d), max(d), len=100)
	w.bgauss <- w_bgauss(ud, L, knots)
	plot(ud, w.bgauss[,1], type="l", ylim=c(0,1));
	for(j in 2:L) { lines(ud,w.bgauss[,j],lty=j) }
}

"w_sgauss" <- function(d, L, knots) {
	# construct weights from Gaussian kernels with varying bandwidths
	if (L > 1) {
		sigmas <- rep(0, length(knots))
		for (i in 1:length(knots)) {
			if (i == 1) {
				sigmas[i] <- 0.25*(knots[i+1]-knots[i])
			} else if (i == length(knots)) {
				sigmas[i] <- 0.25*(knots[i]-knots[i-1])
			} else {
				sigmas[i] <- 0.25*min( c(knots[i]-knots[i-1], knots[i+1]-knots[i]) )
			}
		}

		diffs <- sapply(d, function(i) { i-knots })

		weights <- t(sapply(1:nrow(diffs), function(row) {
			exp(-diffs[row,]^2/(2*sigmas[row]^2))
		}))

		weights <- apply(weights, 2, function(col) { col/sum(col) })
	}

	t(weights)
}

if (FALSE) {
	L <- 4
	#knots <- seq(min(d), max(d), len=L)
	knots <- c(min(d), 0.1, max(d)/2, max(d))
	ud <- seq(min(d), max(d), len=100)
	w.sgauss <- w_sgauss(ud, L, knots)
	plot(ud, w.sgauss[,1], type="l", ylim=c(0,1));
	for(j in 2:L) { lines(ud,w.sgauss[,j],lty=j) }
}

"w_linear" <- function(d, L, knots) {
	# construct weights from linear functions between knots
	if (L > 1) {
		weights <- matrix(0, nrow=length(d), ncol=L)
		for (l in 1:L) {
			if (l == 1) {
				weights[,l] <- as.integer(d >= knots[l]&d <= knots[l+1])*(knots[l+1]-d)/(knots[l+1]-knots[l])
			} else if (l == L) {
				weights[,l] <- as.integer(d >= knots[l-1]&d <= knots[l])*(1-(knots[l]-d)/(knots[l]-knots[l-1]))
			} else {
				weights[,l] <- as.integer(d >= knots[l-1]&d < knots[l])*(1-(knots[l]-d)/(knots[l]-knots[l-1])) +
				               as.integer(d >= knots[l]&d < knots[l+1])*(knots[l+1]-d)/(knots[l+1]-knots[l])
			}
		}
	}

	weights
}

if (TRUE) {
	L <- 4
	#knots <- seq(min(d), max(d), len=L)
	knots <- c(min(d), 0.1, max(d)/2, max(d))
	ud <- unique(as.vector(unlist( sapply(1:(length(knots)-1), function(i) { seq(knots[i],knots[i+1],len=10) }) )))
	w.linear <- w_linear(ud, L, knots)
	plot(ud, w.linear[,1], type="l", ylim=c(0,1));
	for(j in 2:L) { lines(ud,w.linear[,j],lty=j) }
}
