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
	matrix[k,k]                  cholOmega[L];
}
model {
	matrix[n,k] mu;
	matrix[n,k] sigma;

	for (i in 2:n) {
		for (j in 1:k) {
			// fill in mu and sigma
			mu[i,j] <- 0;

			if ( (j+1) < k ) {
				for (l in (j+1):k) {
					mu[i,j] <- mu[i,j] + cholOmega[k][Zstar[i,l]
				}
			}
		}
	}

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
#	matrix[k,k] Omega[L];
#	matrix[k,k] Sigma_d[Nt];
	matrix[k,k] corrSigma_d[Nt];

#	// construct Omega
#	for (l in 1:L)
#		for (i in 1:k)
#			for (j in 1:k)
#				Omega[l,i,j] <- omega[l,i]*omega[l,j]*corrOmega[l,i,j];
#
#	// construct Sigma(T) = sum_l w_l(T) Omega_l
#	for (t in 1:Nt) {
#		Sigma_d[t] <- weights[1,t] * Omega[1];
#
#		for (l in 2:L)
#			Sigma_d[t] <- Sigma_d[t] + weights[l,t] * Omega[l];
#	}
#
#	// construct corrSigma(T)
#	for (t in 1:Nt)
#		for (i in 1:k)
#			for (j in 1:k)
#					corrSigma_d[t,i,j] <- Sigma_d[t,i,j]/sqrt(Sigma_d[t,i,i]*Sigma_d[t,j,j]);

	// construct Sigma(T) = sum_l w_l(T) Omega_l
	for (t in 1:Nt) {
		corrSigma_d[t] <- weights[1,t] * corrOmega[1];

		for (l in 2:L)
			corrSigma_d[t] <- corrSigma_d[t] + weights[l,t] * corrOmega[l];
	}

}
