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
	vector<lower=-2,upper=2>[L]     s[k];
	vector<lower=-2,upper=2>[L]     r[krho];
}

model {
	matrix[k,k] Sigma;
	int c;

	for (i in 1:k)    s[i] ~ cauchy(0, 5);
	for (i in 1:krho) r[i] ~ cauchy(0, 5);
print(s);
print(r);

	// model
	for (i in 2:n) {
		// compute Sigma(h_i)

		for (k1 in 1:k) Sigma[k1,k1] <- pow(exp(dot_product(weights[i], s[k1])), 2);

		c <- 1;
		for (k1 in 1:(k-1)) {
			for (k2 in (k1+1):k) {
				Sigma[k1,k2] <- (2*inv_logit( dot_product(weights[i], r[c]) )-1)*sqrt(Sigma[k1,k1]*Sigma[k2,k2]);
				Sigma[k2,k1] <- Sigma[k1,k2];
				c <- c+1;
			}
		}

print(Sigma);
		Zstar[i] ~ multi_normal(mu_zero, Sigma);
	}
}
