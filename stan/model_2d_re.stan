data {
	int<lower=0>  n;         // sample size
	int<lower=0>  k;         // number of data sources
	vector[k]     Zstar[n];  // normalized data Zstar

	int<lower=0>  L;         // number of cov matrices in mixture

	vector[L]     weights[n];   // weight to each cov mat
	int           i_max_w[n];   // index of max weight

	int<lower=0>  Nud;        // how many unique d_i do we have?
	vector[L]     udw[Nud];   // weights for the unique Zstar
}
transformed data {
	vector[L]     rt_weights[n];
	vector[n]     max_w;        // value of max weight

	for (i in 1:n)
		for (l in 1:L)
			rt_weights[i,l] <- sqrt(weights[i,l]);

	for (i in 1:n)
		max_w[i] <- weights[i,i_max_w[i]];

	print("n = ", n, ", k = ", k, ", L = ", L, ", Nud = ", Nud);
}

parameters {
	vector<lower=0>[k]         omega[L];
	corr_matrix[k]             corrOmega[L];
	vector[k]                  alpha[n,L];
}

transformed parameters {
	matrix[k,k]  Omega[L];
	vector[k]    mu[n];

	// mean for the Zstar's
	for (i in 1:n)
		for (j in 1:k)
			mu[i,j] <- 0.0;

	for (l in 1:L)
		for (i in 1:k)
			for (j in 1:k)
				Omega[l,i,j] <- omega[l,i]*omega[l,j]*corrOmega[l,i,j];

	{
		matrix[k,k]  cholOmega;

		if (L > 1) {
			for (l in 1:L) {
				cholOmega <- cholesky_decompose(Omega[l]);

				for (i in 2:n) {
					if (i_max_w[i] != l) {
						mu[i] <- mu[i] + rt_weights[i,l] * cholOmega * alpha[i,l];
					}
				}
			}
		}

	}
}

model {
	// priors
	for (i in 2:n)
		for (l in 1:L)
			if (i_max_w[i] != l) alpha[i,l] ~ normal(0, 1);

	for (l in 1:L) omega[l] ~ cauchy(0, 5);
	// implied uniform prior on Omega's

	// model
	for (i in 2:n)
		Zstar[i] ~ multi_normal(mu[i], max_w[i] * Omega[i_max_w[i]]);
}

generated quantities {
	matrix[k,k] corrSigma_d[Nud];

	{
		matrix[k,k] wSigma;

		// construct Sigma(d) = sum_l w_l(d) Omega_l
		for (i in 1:Nud) {
			wSigma <- udw[i,1] * Omega[1];

			for (l in 2:L)
				wSigma <- wSigma + udw[i,l] * Omega[l];

			for (j1 in 1:k) {
				for (j2 in 1:k) {
					corrSigma_d[i,j1,j2] <- wSigma[j1,j2] / sqrt( wSigma[j1,j1] * wSigma[j2,j2] );
				}
			}

		}
	}
}
