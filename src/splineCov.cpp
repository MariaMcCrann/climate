#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "model.h"
#include "hmc.h"
#include "inline.h"

#include <R.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>

class ModelSplineCov : public Model {
public:
	class Prior {
	public:
		double sd;   // prior SD of spline coefficients
		double var;
	};

	class Data {
	public:
		int n;                       // number of observations
		int k;                       // number of sources
		const double *y;             // [n,k] observations
		int L;                       // number of basis functions
		const int    *Nnz;           // [n] number of non-zero weights for this obs
		const int    *Mnz;           // [n,4] matrix of the non-zero indices
		const double *Wnz;           // [n,4] matrix of the non-zero weights
	};

	ModelSplineCov(const Prior *prior, const Data *data);
	~ModelSplineCov();

  virtual int  num_params() const;
	virtual bool log_kernel(const double *theta, double *lp, bool do_lik=true, bool do_pri=true);
	virtual bool grad_lk(const double *theta, double *grad);
	virtual bool scales(const double *theta, double *s);

	void fillR_i(int irow, const double *diag, const double *offd);

private:
	const Prior *mPrior;
	const Data  *mData;
	int          mKTri;       // number of upper triangular elements
	int          mNumPerL;    // number of params for each basis function
	int          mNumParams;  // total number of params
	double      *mR_i;        // [k,k] covariance
	double      *mInvR_i;     // [k,k] inverse covariance
};

ModelSplineCov::ModelSplineCov(const Prior *prior, const Data *data) {
	mPrior = prior;
	mData = data;

	mKTri = mData->k*(mData->k-1)/2;
	mNumPerL = mData->k+mKTri;
	mNumParams = mData->L*mNumPerL;

	mR_i = (double *)malloc(sizeof(double)*mData->k*mData->k);
	for (int i = 0; i < mData->k*mData->k; i++) mR_i[i] = 0;

	mInvR_i = (double *)malloc(sizeof(double)*mData->k*mData->k);
	for (int i = 0; i < mData->k*mData->k; i++) mInvR_i[i] = 0;
}

ModelSplineCov::~ModelSplineCov() {
	free(mR_i);
	free(mInvR_i);
}

int ModelSplineCov::num_params() const {
	return(mNumParams);
}

void ModelSplineCov::fillR_i(int irow, const double *diag, const double *offd) {
	int k1,k2;
	int c,i,index;

	// diagonal elements
	for (k1 = 0; k1 < mData->k; k1++) {
		index = k1 + k1*mData->k;
		mR_i[index] = 0;
		for (i = 0; i < mData->Nnz[irow]; i++) {
			mR_i[index] += diag[k1 + mData->Mnz[irow + i*mData->n]*mData->k]*mData->Wnz[irow + i*mData->n];
		}
		mR_i[index] = exp(mR_i[index]);
	}

	// upper triangular elements
	c = 0;
	for (k2 = 1; k2 < mData->k; k2++) {
		for (k1 = 0; k1 < k2; k1++) {
			index = k1 + k2*mData->k;
			mR_i[index] = 0;
			for (i = 0; i < mData->Nnz[irow]; i++) {
				mR_i[index] += offd[c + mData->Mnz[irow + i*mData->n]*mKTri]*mData->Wnz[irow + i*mData->n];
			}

			c++;
		}
	}

}

bool ModelSplineCov::log_kernel(const double *theta, double *lp, bool do_lik, bool do_pri) {
	double llik = 0;
	double lpri = 0;

	int irow;
	int i,j;

	const double *diag = theta;
	const double *offd = theta + mData->k*mData->L;

	char   cside = 'L';
	char   cuplo = 'U';
	char   ctran = 'T';
	char   cdiag = 'N';
	int    ione   = 1;
	double done   = 1;
	double bs[mData->k];

//for (i = 0; i < mData->k*mData->L; i++) MSG("%.2f; ", diag[i]); MSG("\n");
//for (i = 0; i < mKTri*mData->L; i++) MSG("%.2f; ", offd[i]); MSG("\n");

	if (do_lik) { // likelihood
		// compute column sums of diag
		double csd[mData->L];
		for (i = 0; i < mData->L; i++) {
			csd[i] = 0;
			for (j = 0; j < mData->k; j++) {
				csd[i] += diag[j + i*mData->k];
			}
		}

		for (irow = 0; irow < mData->n; irow++) {
			// fill in R_i
			fillR_i(irow, diag, offd);

/*
if (irow == 2) {
MSG("R_i:\n");
for (k1 = 0; k1 < mData->k; k1++) {
	for (k2 = 0; k2 < mData->k; k2++) {
		MSG("%.2f; ", mR_i[k1 + k2*mData->k]);
	}
	MSG("\n");
} }
*/

			// add in -t(weights) colSums(diag) component
			for (i = 0; i < mData->Nnz[irow]; i++) {
				llik -= mData->Wnz[irow + i*mData->n] * csd[mData->Mnz[irow + i*mData->n]];
			}

			// add in -0.5*t(y[i,]) %*% chol2inv(R_i) %*% y[i,] component

			// solve t(R_i) b = y[i,]
			for (i = 0; i < mData->k; i++) { bs[i] = mData->y[irow + i*mData->n]; };
			dtrsm_(&cside, &cuplo, &ctran, &cdiag,
			       /*m*/&mData->k, /*n*/&ione, /*alpha*/&done, /*A*/mR_i,
			       /*lda*/&mData->k, /*B*/bs, /*ldb*/&mData->k);

			for (i = 0; i < mData->k; i++) {
				llik -= 0.5*pow(bs[i],2);
			}
		}
	}


	if (do_pri) { // prior
		for (i = 0; i < mNumParams; i++) {
			lpri += pow(theta[i],2);
		}
		lpri *= -0.5/pow(mPrior->sd, 2);
	}

	// add likelihood and prior components
	*lp = llik + lpri;

	return(true);
}

bool ModelSplineCov::grad_lk(const double *theta, double *grad) {
	int i,j,k;
	int c;
	int iL,irow;
	int k1,k2;
	int index;
	double w;

	const double *diag = theta;
	const double *offd = theta + mData->k*mData->L;

	double *grad_d = grad;
	double *grad_o = grad + mData->k*mData->L;

	char   cside = 'L';
	char   cuplo = 'U';
	char   ctran = 'N';
	char   cdiag = 'N';
	int    ione  = 1;
	double done  = 1;
	double dzero = 0;

	double A[mData->k*mData->k];
	double A_y_i[mData->k];
	double d;
	double P[mData->k*mData->k];

	// start with prior component
	for (i = 0; i < mNumParams; i++) {
		grad[i] = -theta[i]/mPrior->var;
	}

	// add components from likelihood
	for (irow = 0; irow < mData->n; irow++) {

		// fill in R_i
		fillR_i(irow, diag, offd);

		// invert R_i by solving R_i b = diag(k)
		for (i = 0; i < mData->k*mData->k; i++) { mInvR_i[i] = 0; }
		for (i = 0; i < mData->k; i++) { mInvR_i[i + i*mData->k] = 1; }
		cside = 'L'; cuplo = 'U'; ctran = 'N'; cdiag = 'N';
		dtrsm_(&cside, &cuplo, &ctran, &cdiag,
		       /*m*/&mData->k, /*n*/&mData->k, /*alpha*/&done, /*A*/mR_i,
		       /*lda*/&mData->k, /*B*/mInvR_i, /*ldb*/&mData->k);

		// compute A_y_i = invR_i x t(invR_i) x y[i,]
			// multiply A = invR_i x t(invR_i)
			for (i = 0; i < mData->k*mData->k; i++) { A[i] = mInvR_i[i]; }
			cside = 'R'; cuplo = 'U'; ctran = 'T'; cdiag = 'N';
			dtrmm_(&cside, &cuplo, &ctran, &cdiag,
			       /*m*/&mData->k, /*n*/&mData->k, /*alpha*/&done, /*A*/mInvR_i,
			       /*lda*/&mData->k, /*B*/A, /*ldb*/&mData->k);

			// multiple A_y_i = A x y_i
			ctran = 'N';
			dgemv_(&ctran, /*m*/&mData->k, /*n*/&mData->k, /*alpha*/&done, /*A*/A,
			       /*lda*/&mData->k, /*X*/mData->y+irow, /*incx*/&mData->n,
			       /*beta*/&dzero, /*y*/A_y_i, /*incy*/&ione);

		// deterimenant piece for diagonal elements
		for (i = 0; i < mData->Nnz[irow]; i++) {
			index = irow + i*mData->n;
			w = mData->Wnz[index];
			for (k1 = 0; k1 < mData->k; k1++) {
				grad_d[k1 + mData->Mnz[index]*mData->k] -= w;
			}
		}

		// quadratic form piece
		for (i = 0; i < mData->Nnz[irow]; i++) {
			index = irow + i*mData->n;
			iL = mData->Mnz[index];

			// diagonal elements
			for (k1 = 0; k1 < mData->k; k1++) {
				// compute R_i[k1,k1]*weights[irow,iL]
				d = mR_i[k1 + k1*mData->k] * mData->Wnz[index];

				// compute D x R_i, where D is 0 except for D[k1,k1] = d
				for (j = 0; j < mData->k*mData->k; j++) { P[j] = 0; }
				for (j = k1; j < mData->k; j++) {
					P[k1 + j*mData->k] = d * mR_i[k1 + j*mData->k];
				}

				// compute the rest of the quadratic form
				for (k = k1; k < mData->k; k++) {
					d = 0;
					for (j = k1; j < mData->k; j++) {
						d += (P[j + k*mData->k] + P[k + j*mData->k]) * A_y_i[j];
					}
					grad_d[k1 + iL*mData->k] += A_y_i[k]*d/2;
				}
			}

			// upper triangular elements
			c = 0;
			for (k2 = 1; k2 < mData->k; k2++) {
				for (k1 = 0; k1 < k2; k1++) {
					d = mData->Wnz[index];

					// compute D x R_i, where D is 0 except for D[k2,k1] = d
					for (j = 0; j < mData->k*mData->k; j++) { P[j] = 0; }
					for (j = k1; j < mData->k; j++) {
						P[k2 + j*mData->k] = d * mR_i[k1 + j*mData->k];
					}
/*
if (irow+1==12) {
for (int a = 0; a < mData->k; a++) { for (int b = 0; b < mData->k; b++) { MSG("%.2f ", P[a + b*mData->k]); } MSG("\n"); }
}
*/

					// compute the rest of the quadratic form
					for (k = k1; k < mData->k; k++) {
						d = 0;
						for (j = k1; j < mData->k; j++) {
							d += (P[j + k*mData->k] + P[k + j*mData->k]) * A_y_i[j];
						}
						grad_o[c + iL*mKTri] += A_y_i[k]*d/2;
					}

					c++;
				}
			}

		}  // end non-zero elements

	}

	return(true);
}

bool ModelSplineCov::scales(const double *theta, double *s) {
	int i,j,k;
	int c;
	int iL,irow;
	int k1,k2;
	int index;
	double d;

	const double *diag = theta;
	const double *offd = theta + mData->k*mData->L;

	double *s_d = s;
	double *s_o = s + mData->k*mData->L;

	char   cside = 'L';
	char   cuplo = 'U';
	char   ctran = 'N';
	char   cdiag = 'N';
	double done  = 1;
	double dzero = 0;

	double invSigma[mData->k*mData->k];
	double P[mData->k*mData->k];
	double A[mData->k*mData->k];

	// initialize scales
	for (i = 0; i < mNumParams; i++) {
		s[i] = 0;
	}

	// add components from likelihood
	for (irow = 0; irow < mData->n; irow++) {

		// fill in R_i
		fillR_i(irow, diag, offd);

		// invert R_i by solving R_i b = diag(k)
		for (i = 0; i < mData->k*mData->k; i++) { mInvR_i[i] = 0; }
		for (i = 0; i < mData->k; i++) { mInvR_i[i + i*mData->k] = 1; }
		cside = 'L'; cuplo = 'U'; ctran = 'N'; cdiag = 'N';
		dtrsm_(&cside, &cuplo, &ctran, &cdiag,
		       /*m*/&mData->k, /*n*/&mData->k, /*alpha*/&done, /*A*/mR_i,
		       /*lda*/&mData->k, /*B*/mInvR_i, /*ldb*/&mData->k);

		// compute invSigma = invR_i x t(invR_i)
		for (i = 0; i < mData->k*mData->k; i++) { invSigma[i] = mInvR_i[i]; }
		cside = 'R'; cuplo = 'U'; ctran = 'T'; cdiag = 'N';
		dtrmm_(&cside, &cuplo, &ctran, &cdiag,
		       /*m*/&mData->k, /*n*/&mData->k, /*alpha*/&done, /*A*/mInvR_i,
		       /*lda*/&mData->k, /*B*/invSigma, /*ldb*/&mData->k);

		for (i = 0; i < mData->Nnz[irow]; i++) {
			index = irow + i*mData->n;
			iL = mData->Mnz[index];

			// diagonal elements
			for (k1 = 0; k1 < mData->k; k1++) {
				// compute R_i[k1,k1]*weights[irow,iL]
				d = mR_i[k1 + k1*mData->k] * mData->Wnz[index];

				// compute D x R_i, where D is 0 except for D[k1,k1] = d
				for (j = 0; j < mData->k*mData->k; j++) { P[j] = 0; }
				for (j = k1; j < mData->k; j++) {
					P[k1 + j*mData->k] = d * mR_i[k1 + j*mData->k];
					P[j + k1*mData->k] += P[k1 + j*mData->k];
				}

				// compute A = invSigma x P
				ctran = 'N';
				dgemm_(&ctran, &ctran,
				       /*m*/&mData->k, /*n*/&mData->k, /*k*/&mData->k, /*alpha*/&done,
				       /*A*/invSigma, /*lda*/&mData->k, /*B*/P, /*ldb*/&mData->k,
				       /*beta*/&dzero, /*C*/A, /*ldc*/&mData->k);

				// add trace to scale
				for (k = 0; k < mData->k; k++) {
					for (j = 0; j < mData->k; j++) {
						s_d[k1 + iL*mData->k] += A[k + j*mData->k] * A[j + k*mData->k];
					}
				}

			}

			// upper triangular elements
			c = 0;
			for (k2 = 1; k2 < mData->k; k2++) {
				for (k1 = 0; k1 < k2; k1++) {
					d = mData->Wnz[index];

					// compute t(D) x R_i, where D is 0 except for D[k1,k2] = d
					for (j = 0; j < mData->k*mData->k; j++) { P[j] = 0; }
					for (j = 0; j < mData->k; j++) {
						P[k2 + j*mData->k] = d * mR_i[k1 + j*mData->k];
						P[j + k2*mData->k] += P[k2 + j*mData->k];
					}

					// compute A = invSigma x P
					ctran = 'N';
					dgemm_(&ctran, &ctran,
					       /*m*/&mData->k, /*n*/&mData->k, /*k*/&mData->k, /*alpha*/&done,
					       /*A*/invSigma, /*lda*/&mData->k, /*B*/P, /*ldb*/&mData->k,
					       /*beta*/&dzero, /*C*/A, /*ldc*/&mData->k);

					// add trace to scale
					for (k = 0; k < mData->k; k++) {
						for (j = 0; j < mData->k; j++) {
							s_o[c + iL*mKTri] += A[k + j*mData->k] * A[j + k*mData->k];
						}
					}

					c++;
				}
			}

		}
	}

	// SD scale
	for (i = 0; i < mNumParams; i++) {
		s[i] = 1/sqrt(0.5*s[i] + 1/mPrior->var);
	}

	return(true);
}

extern "C" {

void spline_cov_fit(
	// prior
	double *prior,
	// data
	int *n, int *k,
	double *y,
	int *L,
	int *Nnz, int *Mnz, double *Wnz,
	// sampler 
	double *step_e, int *step_L,
	double *inits,
	int *Niter, double *samples,
	bool *verbose
) {

	// prior
	ModelSplineCov::Prior *model_prior = new ModelSplineCov::Prior();
	model_prior->sd  = prior[0];
	model_prior->var = pow(prior[0], 2);

	// data
	ModelSplineCov::Data *model_data = new ModelSplineCov::Data();
	model_data->n   = n[0];
	model_data->k   = k[0];
	model_data->y   = (const double *)y;
	model_data->L   = L[0];
	model_data->Nnz = (const int *)Nnz;
	model_data->Mnz = (const int *)Mnz;
	model_data->Wnz = (const double *)Wnz;

	// model
	Model *m = new ModelSplineCov((const ModelSplineCov::Prior *)model_prior, (const ModelSplineCov::Data *)model_data);

/*
	double lp;
	m->log_kernel((const double *)inits, &lp);
	//MSG("lp = %.2f\n", lp);

	double grad[m->num_params()];
	m->grad_lk((const double*)inits, grad);
	//for (int i = 0; i < m->num_params(); i++) MSG("%.2f ", grad[i]); MSG("\n");

	double scales[m->num_params()];
	m->scales((const double*)inits, scales);
	for (int i = 0; i < m->num_params(); i++) MSG("%.2f ", scales[i]); MSG("\n");
*/

	// get samples
	HMC *sampler = new HMC(m, (const double *)inits, samples);
	sampler->sample(*Niter, *step_e, *step_L, verbose);

	delete sampler;
	delete m;
	delete model_prior;
	delete model_data;
}

void spline_cov_lk(
	// prior
	double *prior,
	// data
	int *n, int *k,
	double *y,
	int *L,
	int *Nnz, int *Mnz, double *Wnz,
	// evaluation point
	double *eval,
	// result
	double *lk
) {

	// prior
	ModelSplineCov::Prior *model_prior = new ModelSplineCov::Prior();
	model_prior->sd  = prior[0];
	model_prior->var = pow(prior[0], 2);

	// data
	ModelSplineCov::Data *model_data = new ModelSplineCov::Data();
	model_data->n   = n[0];
	model_data->k   = k[0];
	model_data->y   = (const double *)y;
	model_data->L   = L[0];
	model_data->Nnz = (const int *)Nnz;
	model_data->Mnz = (const int *)Mnz;
	model_data->Wnz = (const double *)Wnz;

	// model
	Model *m = new ModelSplineCov((const ModelSplineCov::Prior *)model_prior, (const ModelSplineCov::Data *)model_data);

	m->log_kernel((const double *)eval, lk);

	delete m;
	delete model_prior;
	delete model_data;
}

void spline_cov_gr(
	// prior
	double *prior,
	// data
	int *n, int *k,
	double *y,
	int *L,
	int *Nnz, int *Mnz, double *Wnz,
	// evaluation point
	double *eval,
	// result
	double *gr
) {

	// prior
	ModelSplineCov::Prior *model_prior = new ModelSplineCov::Prior();
	model_prior->sd  = prior[0];
	model_prior->var = pow(prior[0], 2);

	// data
	ModelSplineCov::Data *model_data = new ModelSplineCov::Data();
	model_data->n   = n[0];
	model_data->k   = k[0];
	model_data->y   = (const double *)y;
	model_data->L   = L[0];
	model_data->Nnz = (const int *)Nnz;
	model_data->Mnz = (const int *)Mnz;
	model_data->Wnz = (const double *)Wnz;

	// model
	Model *m = new ModelSplineCov((const ModelSplineCov::Prior *)model_prior, (const ModelSplineCov::Data *)model_data);

	m->grad_lk((const double*)eval, gr);

	delete m;
	delete model_prior;
	delete model_data;
}

} // end extern "C"
