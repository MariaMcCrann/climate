// unconstrained/Cholesky factorization
data {
	int<lower=0>  n;           // sample size
	int<lower=0>  k;           // number of data sources
	vector[k]     Zstar[n];    // normalized data Zstar

	int<lower=0>  L;           // number of variables in mixture

	vector[L]     weights[n];  // weight to each cov mat
	int           i_max_w[n];  // index of max weight

	int   nz_max;              // max number non-zero
	int   Nnz[n];              // number of non-zero elements for this observation
	int   Mnz[n,nz_max];       // matrix with non-zero indices
	vector[nz_max] Wnz[n];     // matrix with non-zero weights

	int<lower=0>  Nuf;         // how many unique f_i do we have?
	vector[L]     ufw[Nuf];    // weights for the unique Zstar

	int krho;
}

transformed data {
	vector[k]     mu_zero;

	for (i in 1:k) mu_zero[i] <- 0.0;

	print("n = ", n, ", k = ", k, ", L = ", L, ", krho = ", krho);
}

parameters {
	vector[L]     s[k];
	vector[L]     r[krho];
}

model {
	matrix[k,k] SigmaL;
	int c;

	// priors
	for (i in 1:k)    s[i] ~ cauchy(0, 5);
	for (i in 1:krho) r[i] ~ cauchy(0, 5);

	// model
	for (i in 2:n) {
		// compute Sigma(h_i)

		for (k1 in 1:k) {
			SigmaL[k1,k1] <- 0;
			for (j in 1:Nnz[i]) {
				SigmaL[k1,k1] <- SigmaL[k1,k1] + Wnz[i,j] * s[k1,Mnz[i,j]];
			}
			SigmaL[k1,k1] <- exp(SigmaL[k1,k1]);
		}

		c <- 1;
		for (k1 in 1:(k-1)) {
			for (k2 in (k1+1):k) {
				SigmaL[k1,k2] <- 0;
				SigmaL[k2,k1] <- 0;

				for (j in 1:Nnz[i]) {
					SigmaL[k2,k1] <- SigmaL[k2,k1] + Wnz[i,j] * r[c,Mnz[i,j]];
				}
				c <- c+1;
			}
		}

		Zstar[i] ~ multi_normal_cholesky(mu_zero, SigmaL);
	}
}
