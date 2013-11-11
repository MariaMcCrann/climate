data {
	int<lower=0>  n;         // sample size
	int<lower=0>  k;         // number of data sources
	vector[k]     Zstar[n];  // normalized data Zstar

	int<lower=0>  L;         // number of cov matrices in mixture

	vector[L]     weights[n];   // weight to each cov mat
	int           i_max_w[n];   // index of max weight

	int<lower=0>  Nuf;        // how many unique d_i do we have?
	vector[L]     ufw[Nuf];   // weights for the unique Zstar
}

transformed data {
	vector[L]     rt_weights[n];
	vector[n]     max_w;        // value of max weight
	vector[k]     mu;

	// mean for the Zstar's
	for (i in 1:k)
		mu[i] <- 0.0;

	for (i in 1:n)
		for (l in 1:L)
			rt_weights[i,l] <- sqrt(weights[i,l]);

	for (i in 1:n)
		max_w[i] <- weights[i,i_max_w[i]];

	print("n = ", n, ", k = ", k, ", L = ", L, ", Nuf = ", Nuf);
}

parameters {
	vector<lower=0>[k]         omega[L];
	corr_matrix[k]             corrOmega[L];
}

transformed parameters {
	matrix[k,k]  Omega[L];

	for (l in 1:L)
		for (i in 1:k)
			for (j in 1:k)
				Omega[l,i,j] <- omega[l,i]*omega[l,j]*corrOmega[l,i,j];
}

model {
	matrix[k,k] wSigma;

	for (l in 1:L) omega[l] ~ cauchy(0, 5);
	// implied uniform prior on Omega's

	// model
	for (i in 2:n) {
		wSigma <- weights[i,1] * Omega[1];

		for (l in 2:L)
			wSigma <- wSigma + weights[i,l] * Omega[l];

		Zstar[i] ~ multi_normal(mu, wSigma);
	}
}

generated quantities {
	matrix[k,k] corrSigma_f[Nuf];

	{
		matrix[k,k] wSigma;

		// construct Sigma(f) = sum_l w_l(f) Omega_l
		for (i in 1:Nuf) {
			wSigma <- ufw[i,1] * Omega[1];

			for (l in 2:L)
				wSigma <- wSigma + ufw[i,l] * Omega[l];

			for (j1 in 1:k) {
				for (j2 in 1:k) {
					corrSigma_f[i,j1,j2] <- wSigma[j1,j2] / sqrt( wSigma[j1,j1] * wSigma[j2,j2] );
				}
			}

		}
	}
}
