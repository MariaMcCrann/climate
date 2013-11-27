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
}

transformed data {
	vector[L]     rt_weights[n];
	vector[n]     max_w;        // value of max weight
	vector[k]     mu_zero;
	matrix[k,k]   zeroMat;
	real ksq;
	int krho;

	for (i in 1:n)
		for (l in 1:L)
			rt_weights[i,l] <- sqrt(weights[i,l]);

	for (i in 1:n)
		max_w[i] <- weights[i,i_max_w[i]];

	for (i in 1:k) for (j in 1:k) zeroMat[i,j] <- 0.0;

	for (i in 1:k) mu_zero[i] <- 0.0;

	ksq <- k*k;
	krho <- k*(k-1)/2;

	print("n = ", n, ", k = ", k, ", L = ", L, ", Nuf = ", Nuf);
}

parameters {
	vector[L]     v[k];
	vector[L]     r[krho];
}

model {
	matrix[k,k] Sigma;
	vector[k]   sigma2;
	//matrix[k,k] rho;
	int c;

	for (i in 1:k)    v[i] ~ cauchy(0, 5);
	for (i in 1:krho) r[i] ~ cauchy(0, 5);

	// model
	for (i in 2:n) {
		// compute Sigma(h_i)

/*
		for (m in 1:k) {
			sigma2[m] <- 0;
			for (j in 1:Nnz[i]) sigma2[m] <- sigma2[m] + weights[i,Mnz[i,j]] * v[k][Mnz[i,j]];
			sigma2[m] <- exp(sigma2[m]);
		}

		sigma2 <- exp(sigma2);
*/

		c <- 1;
		for (k1 in 1:k) {
			sigma2[k1] <- exp(dot_product(weights[i] , v[k]));
		}

		for (k1 in 1:k) {
			Sigma[k1,k1] <- sigma2[k1];

			for (k2 in (k1+1):k) {
				Sigma[k1,k2] <- (2*inv_logit( dot_product(weights[i], r[c]) )-1)*sqrt(sigma2[k1]*sigma2[k2]);
				Sigma[k2,k1] <- Sigma[k1,k2];
				c <- c+1;
			}
		}

		Zstar[i] ~ multi_normal(mu_zero, Sigma);
	}
}
