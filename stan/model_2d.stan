data {
	int<lower=0>  n;         // sample size
	int<lower=0>  k;         // number of data sources
	vector[k]     Zstar[n];  // normalized data Zstar

	int<lower=0>  L;         // number of cov matrices in mixture

	vector[L]     weights[n];   // weight to each cov mat
	int           i_max_w[n];   // index of max weight
}

transformed data {
	vector[k]     mu_zero;

	for (i in 1:k) mu_zero[i] <- 0;
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
	// local variables
	matrix[k,k] wSigma;

	// priors
	for (l in 1:L) omega[l] ~ cauchy(0, 5);
	// implied uniform prior on Omega's

	// model
	for (i in 2:n) {
		// construct Sigma(f) = sum_l w_l(j) Omega_l
		wSigma <- weights[i,1] * Omega[1];
		for (l in 2:L) wSigma <- wSigma + weights[i,l] * Omega[l];

		Zstar[i] ~ multi_normal(mu_zero, wSigma);
	}
}

generated quantities {
	real Dbar;

	{
		matrix[k,k] wSigma;
		vector[n] dic;

		// compute Dbar
		dic[1] <- 0.0;

		for (i in 2:n) {
			wSigma <- weights[i,1] * Omega[1];
			for (l in 2:L) wSigma <- wSigma + weights[i,l] * Omega[l];

			dic[i] <- multi_normal_log(Zstar[i], mu_zero, wSigma);
		}
		Dbar <- -2*sum(dic);
	}
}
