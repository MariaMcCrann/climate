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

	int<lower=0>  Nuf;         // how many unique f_i do we have?
	vector[L]     ufw[Nuf];    // weights for the unique Zstar

	int krho;
}

transformed data {
	vector[L]     rt_weights[n];
	vector[n]     max_w;        // value of max weight
	vector[k]     mu_zero;
	matrix[k,k]   zeroMat;
	real ksq;

	for (i in 1:n)
		for (l in 1:L)
			rt_weights[i,l] <- sqrt(weights[i,l]);

	for (i in 1:n)
		max_w[i] <- weights[i,i_max_w[i]];

	for (i in 1:k) for (j in 1:k) zeroMat[i,j] <- 0.0;

	for (i in 1:k) mu_zero[i] <- 0.0;

	ksq <- k*k;

	print("n = ", n, ", k = ", k, ", L = ", L, ", krho = ", krho);
}

parameters {
	vector[L]     s[k];
	vector[L]     r[krho];
}

model {
	matrix[k,k] SigmaL;
	int c;

	for (i in 1:k)    s[i] ~ cauchy(0, 5);
	for (i in 1:krho) r[i] ~ cauchy(0, 5);

	// model
	for (i in 2:n) {
		// compute Sigma(h_i)

		for (k1 in 1:k) SigmaL[k1,k1] <- exp(dot_product(weights[i], s[k1]));

		c <- 1;
		for (k1 in 1:(k-1)) {
			for (k2 in (k1+1):k) {
				SigmaL[k2,k1] <- dot_product(weights[i], r[c]);
				SigmaL[k1,k2] <- 0;
				c <- c+1;
			}
		}

		Zstar[i] ~ multi_normal_cholesky(mu_zero, SigmaL);
	}
}

generated quantities {
	vector[k]   sdSigma_f[Nuf];
	matrix[k,k] corrSigma_f[Nuf];
	real        Dbar;

	{
		matrix[k,k] SigmaL;
		matrix[k,k] Sigma;
		vector[n] dic;
		int c;

		// construct variance vector and correlation matrix for Sigma(f_i)
		for (i in 1:Nuf) {
			for (k1 in 1:k) SigmaL[k1,k1] <- exp(dot_product(ufw[i], s[k1]));

			c <- 1;
			for (k1 in 1:(k-1)) {
				for (k2 in (k1+1):k) {
					SigmaL[k2,k1] <- dot_product(ufw[i], r[c]);
					SigmaL[k1,k2] <- 0;
					c <- c+1;
				}
			}

			Sigma <- SigmaL * SigmaL';

			for (k1 in 1:k) {
				sdSigma_f[i,k1] <- sqrt(Sigma[k1,k1]);
				corrSigma_f[i,k1,k1] <- 1;
			}

			for (k1 in 1:(k-1)) {
				for (k2 in (k1+1):k) {
					corrSigma_f[i,k1,k2] <- Sigma[k1,k2] / ( sdSigma_f[i,k1] * sdSigma_f[i,k2] );
					corrSigma_f[i,k2,k1] <- corrSigma_f[i,k1,k2];
				}
			}
		}

		// compute Dbar
		dic[1] <- 0.0;

		for (i in 2:n) {
			for (k1 in 1:k) SigmaL[k1,k1] <- exp(dot_product(weights[i], s[k1]));

			c <- 1;
			for (k1 in 1:(k-1)) {
				for (k2 in (k1+1):k) {
					SigmaL[k2,k1] <- dot_product(weights[i], r[c]);
					SigmaL[k1,k2] <- 0;
					c <- c+1;
				}
			}

			dic[i] <- multi_normal_cholesky_log(Zstar[i], mu_zero, SigmaL);
		}

		Dbar <- -2*sum(dic);
	}
}
