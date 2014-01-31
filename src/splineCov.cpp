#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "model.h"
#include "mh.h"
#include "utils.h"

#include <R.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>

// Specifies covariance over time model:
//
//     y_i   ind N(0, Sigma_i)
// theta_jkl ind N(mu_jk, tau_jk^2)
//    mu_jk  iid N(0, mu_sd^2)
//   tau_jk  iid InvGamma(tau_alpha, tau_beta)
//
// - Sigma_i = R_i'R_i
// - R_i,jj  = exp(sum_l w_i,l theta_jjl)
// - R_i,jk  = sum_l w_i,l theta_jkl        j < k
class ModelSplineCov : public Model {
	friend class SplineCovSampler;
public:
	class Prior {
	public:
		double mu_sd;
		double mu_var;

		double tau_alpha;
		double tau_beta;

// old; not used anymore
		double sd;
		double var;
	};

	class Data {
	public:
		int           n;             // number of observations
		int           k;             // number of sources
		const double *y;             // [n,k] observations
		int           L;             // number of basis functions
		const int    *Nnz;           // [n] number of non-zero weights for this obs
		const int    *Mnz;           // [n,4] matrix of the non-zero indices
		const double *Wnz;           // [n,4] matrix of the non-zero weights
	};

	ModelSplineCov(const Prior *prior, const Data *data);
	~ModelSplineCov();

  virtual int  num_params() const;
	virtual bool log_kernel(const double *theta, double *lp, double *llik=NULL, double *lpri=NULL, bool do_lik=true, bool do_pri=true, int param=-1);
	virtual bool grad_lk(const double *theta, double *grad, const int row=-1);
	virtual bool scales(const double *theta, double *s);
	virtual bool get_obs_info(const double *theta, double *info);

	void fillR_i(int irow, const double *diag, const double *offd);

	int num_theta() { return(mNtheta); }

protected:
	virtual bool obs_info(const double *theta, double *gr);

	const Prior *mPrior;
	const Data  *mData;
	int          mKTri;     // number of upper triangular elements
	int          mNPerL;    // number of params for each basis function
	int          mNtheta;   // number of theta's
	int          mNparams;  // total number of params
	double      *mR_i;      // [k,k] covariance
	double      *mInvR_i;   // [k,k] inverse covariance
	double      *mObsInfo;  // [mNparams,mNparams] observed Fisher information
};

ModelSplineCov::ModelSplineCov(const Prior *prior, const Data *data) {
	mPrior = prior;
	mData = data;

	mKTri    = mData->k*(mData->k-1)/2;
	mNPerL   = mData->k+mKTri;
	mNtheta  = mData->L*mNPerL;
	mNparams = mNtheta + 2*mData->k + 2*mKTri; // counts theta, mu, and tau

	mR_i = (double *)malloc(sizeof(double)*mData->k*mData->k);
	for (int i = 0; i < mData->k*mData->k; i++) mR_i[i] = 0;

	mInvR_i = (double *)malloc(sizeof(double)*mData->k*mData->k);
	for (int i = 0; i < mData->k*mData->k; i++) mInvR_i[i] = 0;

	mObsInfo = (double *)malloc(sizeof(double)*mNparams*mNparams);
	for (int i = 0; i < mNparams*mNparams; i++) mObsInfo[i] = 0;
}

ModelSplineCov::~ModelSplineCov() {
	free(mR_i);
	free(mInvR_i);
	free(mObsInfo);
}

int ModelSplineCov::num_params() const {
	return(mNparams);
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

bool ModelSplineCov::get_obs_info(const double *theta, double *info) {
	int i,j;
	double grad[mNparams*mNparams];

	// get fisher info
	obs_info(theta, grad);

	// add prior piece for diagonals
	for (i = 0; i < mNparams; i++) {
		mObsInfo[i + i*mNparams] += 1/mPrior->var;
	}

	// invert info
	if (chol2inv(mNparams, mObsInfo) != 0) return(false);

	// factorize covariance
	if (chol(mNparams, mObsInfo) != 0) return(false);

	for (i = 0; i < mNparams; i++) {
		for (j = i; j < mNparams; j++) {
			info[i + j*mNparams] = mObsInfo[i + j*mNparams];
		}
	}

	return(true);
}

bool ModelSplineCov::obs_info(const double *theta, double *gr) {
	int irow;
	int i,j;

	// initialize observed info
	for (i = 0; i < mNparams*mNparams; i++) mObsInfo[i] = 0;

	// initialize gradient
	for (i = 0; i < mNparams; i++) gr[i] = 0;

	double grad[mNparams];
	for (irow = 0; irow < mData->n; irow++) {
		grad_lk(theta, grad, irow);

		// add to info
		for (i = 0; i < mNparams; i++) {
			gr[i] += grad[i];

			for (j = i; j < mNparams; j++) {
				mObsInfo[i + j*mNparams] += grad[i]*grad[j];
			}
		}
	}

	return(true);
}

bool ModelSplineCov::log_kernel(const double *theta, double *lk, double *llik, double *lpri, bool do_lik, bool do_pri, int param) {
	double lik = 0;
	double pri = 0;

	int irow;
	int i,j;
	bool do_comp;
	int param_col;

	const double *diag  = theta;
	const double *offd  = diag  + mData->k*mData->L;
	const double *mu_d  = offd  + mKTri*mData->L;
	const double *mu_o  = mu_d  + mData->k;
	const double *tau_d = mu_o  + mKTri;
	const double *tau_o = tau_d + mData->k;

	char   cside = 'L';
	char   cuplo = 'U';
	char   ctran = 'T';
	char   cdiag = 'N';
	int    ione   = 1;
	double done   = 1;
	double bs[mData->k];

	if (param >= 0) {
		// which weight column does this belong to?
		if (param < mData->k*mData->L) { // diagonal
			param_col = floor((double)param/mData->k);
		} else { // upper triangular
			param_col = floor((double)(param-mData->k*mData->L)/mKTri);
		}
//MSG("param = %d, param_col = %d\n", param, param_col);
	} else {
		param_col = 0; // unused
	}

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
			do_comp = true;

			if (param >= 0) {
				// we only care about a single parameter
				do_comp = false;
				for (i = 0; i < mData->Nnz[irow]; i++) {
					if (param_col == mData->Mnz[irow + i*mData->n]) {
						// this weight is non-zero, so add it in
						do_comp = true;
						break;
					}
				}
			}

			if (!do_comp) continue; // skip this row

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
				lik -= mData->Wnz[irow + i*mData->n] * csd[mData->Mnz[irow + i*mData->n]];
			}

			// add in -0.5*t(y[i,]) %*% chol2inv(R_i) %*% y[i,] component

			// solve t(R_i) b = y[i,]
			for (i = 0; i < mData->k; i++) { bs[i] = mData->y[irow + i*mData->n]; };
			dtrsm_(&cside, &cuplo, &ctran, &cdiag,
			       /*m*/&mData->k, /*n*/&ione, /*alpha*/&done, /*A*/mR_i,
			       /*lda*/&mData->k, /*B*/bs, /*ldb*/&mData->k);

			for (i = 0; i < mData->k; i++) {
				lik -= 0.5*pow(bs[i],2);
			}
		}
	}


	if (do_pri) { // prior
		// prior is: p(theta, mu, tau)
		if (param >= 0) {
MSG("Complete prior for single param\n");
			// only add for this parameter
			//pri += pow(theta[param], 2);
			//pri *= -0.5/mPrior->var;
		} else {
			// add for theta
			for (j = 0; j < mData->L; j++) {
				for (i = 0; i < mData->k; i++) pri += -0.5*pow(diag[i + j*mData->k]-mu_d[i],2)/tau_d[i];
				for (i = 0; i < mKTri;    i++) pri += -0.5*pow(offd[i + j*mKTri   ]-mu_o[i],2)/tau_o[i];
			}

			// add for mu
			for (i = 0; i < mData->k; i++) pri += -0.5*pow(mu_d[i],2)/mPrior->mu_var;
			for (i = 0; i < mKTri;    i++) pri += -0.5*pow(mu_o[i],2)/mPrior->mu_var;

			// add for tau
			for (i = 0; i < mData->k; i++) pri += -(mPrior->tau_alpha+1)*log(tau_d[i]) -mPrior->tau_beta/tau_d[i];
			for (i = 0; i < mKTri;    i++) pri += -(mPrior->tau_alpha+1)*log(tau_o[i]) -mPrior->tau_beta/tau_o[i];
		}
	}

	// add likelihood and prior components
	*lk = lik + pri;

	if (llik != NULL) *llik = lik;
	if (lpri != NULL) *lpri = pri;

	return(true);
}

bool ModelSplineCov::grad_lk(const double *theta, double *grad, int row) {
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

	int row_l = 0;
	int row_u = mData->n;

	if (row >= 0) {
		// use a specific row only from likelihood
		row_l = row;
		row_u = row+1;
		for (i = 0; i < mNparams; i++) {
			grad[i] = 0;
		}
	} else {
		// start with prior component
		for (i = 0; i < mNparams; i++) {
			grad[i] = -theta[i]/mPrior->var;
		}
	}

	// add components from likelihood
	for (irow = row_l; irow < row_u; irow++) {

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
	int i;
	double grad[mNparams];

	// use observed info for scales
	obs_info(theta, grad);

	// add prior piece for diagonals
	for (i = 0; i < mNparams; i++) {
		mObsInfo[i + i*mNparams] += 1/mPrior->var;
	}

	// invert info
	if (chol2inv(mNparams, mObsInfo) != 0) return(false);

	// set scales
	for (i = 0; i < mNparams; i++) {
		s[i] = sqrt(mObsInfo[i + i*mNparams]);
	}

/*
	// adjust scales for steep gradient
	for (i = 0; i < mNparams; i++) {
		if (s[i] > 1) { s[i] = 1; }

		if (s[i]*abs(grad[i]) > 0.5) {
			s[i] = 0.5/abs(grad[i]);
		}
	}
*/

	return(true);
}

// do the sampling for this model
class SplineCovSampler {
public:
	SplineCovSampler(ModelSplineCov *m);
	~SplineCovSampler();

	//MH *sampler = new MH(m, (const double *)inits, samples, deviance);
	//sampler->sample(*Niter, *thin, *step_e, true, *verbose);
	bool sample(double *samples, double *deviance, int niter, const double *inits, double e, int nburn=0, bool verbose=false);

private:
	ModelSplineCov *mModel;
	double *mPos;
	double *mPosCur;
	int mIter;
	int mNparams;
	int mNtheta;

	MHSingle **mSamplerTheta;
};

SplineCovSampler::SplineCovSampler(ModelSplineCov *m) {
	mModel   = m;
	mNparams = m->num_params();
	mNtheta  = m->num_theta();

	mPos = NULL;
	mPosCur = NULL;
}

SplineCovSampler::~SplineCovSampler() {
	free(mPos);
}

bool SplineCovSampler::sample(double *samples, double *deviance,
                              int niter, const double *inits, double e, int nburn, bool verbose) {
	int iter;
	int i,p;
	double curLK,curLLik;
	bool burn = true;

	double *diag,*diag_cur;
	double *offd,*offd_cur;
	double *mu_d,*mu_d_cur;
	double *mu_o,*mu_o_cur;
	double *tau_d,*tau_d_cur;
	double *tau_o,*tau_o_cur;

	// for updates
	double theta_mu = 0;
	double update_mu = 0;
	double update_sd = 0;

	// allocate space
	free(mPos);
	mPos = (double *)malloc(sizeof(double)*mNparams);
	free(mPosCur);
	mPosCur = (double *)malloc(sizeof(double)*mNparams);

	diag  = mPos;
	offd  = diag  + mModel->mData->k*mModel->mData->L;
	mu_d  = offd  + mModel->mKTri*mModel->mData->L;
	mu_o  = mu_d  + mModel->mData->k;
	tau_d = mu_o  + mModel->mKTri;
	tau_o = tau_d + mModel->mData->k;

	diag_cur  = mPosCur;
	offd_cur  = diag_cur  + mModel->mData->k*mModel->mData->L;
	mu_d_cur  = offd_cur  + mModel->mKTri*mModel->mData->L;
	mu_o_cur  = mu_d_cur  + mModel->mData->k;
	tau_d_cur = mu_o_cur  + mModel->mKTri;
	tau_o_cur = tau_d_cur + mModel->mData->k;

	// create samples for theta
	mSamplerTheta = new MHSingle*[mNtheta];
	for (p = 0; p < mNtheta; p++) {
		mSamplerTheta[p] = new MHSingle(mModel, p, inits[p], e);
	}

	// set inits
	for (i = 0; i < mNparams; i++) mPos[i] = mPosCur[i] = inits[i];

	// get initial lik
	mModel->log_kernel(mPos, &curLK, &curLLik);
//MSG("%.4f, %.4f\n", curLK, curLLik); niter = 0;

	for (iter = 0; iter < niter; iter++) {
		if (nburn < iter) burn=true;
		else              burn=false;

		// update theta
		for (p = 0; p < mNtheta; p++) {
			mSamplerTheta[p]->sample(mPos, mPosCur, &curLK, &curLLik, burn, verbose);
		}

		// update mu
			// ... diagonal
			for (p = 0; p < mModel->mData->k; p++) {
				// mean of theta_jkl
				theta_mu = 0;
				for (i = 0; i < mModel->mData->L; i++) {
					theta_mu += diag_cur[p + i*mModel->mData->k];
				}
				theta_mu /= mModel->mData->L;

				// update prameters
				update_sd = 1/(mModel->mData->L/tau_d_cur[p] + 1/mModel->mPrior->mu_var);
				update_mu = (mModel->mData->L*theta_mu/tau_d_cur[p])*update_sd;
				update_sd = sqrt(update_sd);

				GetRNGstate();
				mu_d_cur[p] = rnorm(update_mu, update_sd);
				PutRNGstate();
			}
			// ... upper triangle
			for (p = 0; p < mModel->mKTri; p++) {
				// mean of theta_jkl
				theta_mu = 0;
				for (i = 0; i < mModel->mData->L; i++) {
					theta_mu += offd_cur[p + i*mModel->mData->k];
				}
				theta_mu /= mModel->mData->L;

				// update prameters
				update_sd = 1/(mModel->mData->L/tau_o_cur[p] + 1/mModel->mPrior->mu_var);
				update_mu = (mModel->mData->L*theta_mu/tau_o_cur[p])*update_sd;
				update_sd = sqrt(update_sd);

				GetRNGstate();
				mu_o_cur[p] = rnorm(update_mu, update_sd);
				PutRNGstate();
			}

//> a<-3;b<-3;b/(a-1);b^2/( (a-1)^2*(a-2) );v<-1/rgamma(500000,shape=a,scale=1/b);mean(v);var(v);
		// update tau
			// ... diagonal
			for (p = 0; p < mModel->mData->k; p++) {
				// sum of squared differences
				theta_mu = 0;
				for (i = 0; i < mModel->mData->L; i++) {
					theta_mu += pow(diag_cur[p + i*mModel->mData->k] - mu_d[p], 2.0);
				}

				// update prameters
				update_mu = (mModel->mPrior->tau_alpha + 0.5*mModel->mData->L);
				update_sd = 1/(mModel->mPrior->tau_beta + 0.5*theta_mu);

				GetRNGstate();
				tau_d_cur[p] = 1/rgamma(update_mu, update_sd);
				PutRNGstate();
			}
			// ... upper triangle
			for (p = 0; p < mModel->mKTri; p++) {
				// sum of squared differences
				theta_mu = 0;
				for (i = 0; i < mModel->mData->L; i++) {
					theta_mu += pow(offd_cur[p + i*mModel->mData->k] - mu_o[p], 2.0);
				}

				// update prameters
				update_mu = (mModel->mPrior->tau_alpha + 0.5*mModel->mData->L);
				update_sd = 1/(mModel->mPrior->tau_beta + 0.5*theta_mu);

				GetRNGstate();
				tau_o_cur[p] = 1/rgamma(update_mu, update_sd);
				PutRNGstate();
			}

		// save samples
		for (i = 0; i < mNparams; i++) samples[iter + i*niter] = mPosCur[i];

		// get deviance
		mModel->log_kernel(mPosCur, &curLK, &curLLik);
		deviance[iter] = -2*curLLik;

		if (verbose && (iter+1) % 1 == 0) {
			MSG("[%d,%.2f]:", iter+1, curLK);
			for (i = 0; i < 5; i++)
				MSG(" %.2f/%.2f/%.4f", mSamplerTheta[p]->get_local_accept(), mSamplerTheta[p]->get_accept(), mSamplerTheta[p]->get_eps());
			MSG("\n");
		}

#ifndef MATHLIB_STANDALONE
		R_CheckUserInterrupt();
#endif
	}

	for (p = 0; p < mNtheta; p++) {
		delete mSamplerTheta[p];
	}
	delete mSamplerTheta;

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
	double *step_e, int *step_L, int *thin,
	double *inits,
	int *Niter, double *samples, double *deviance,
	bool *verbose
) {

	// prior
	ModelSplineCov::Prior *model_prior = new ModelSplineCov::Prior();
	//model_prior->sd  = prior[0];
	//model_prior->var = pow(prior[0], 2);
	model_prior->mu_sd  = prior[0];
	model_prior->mu_var = pow(prior[0], 2.0);
	model_prior->tau_alpha = prior[1];
	model_prior->tau_beta  = prior[2];

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
	ModelSplineCov *m = new ModelSplineCov((const ModelSplineCov::Prior *)model_prior, (const ModelSplineCov::Data *)model_data);

/*
	double lp;
	m->log_kernel((const double *)inits, &lp);
	//MSG("lp = %.2f\n", lp);

	double grad[m->num_params()];
	m->grad_lk((const double*)inits, grad);
	for (int i = 0; i < m->num_params(); i++) MSG("%.2f ", grad[i]); MSG("\n");

	double scales[m->num_params()];
	m->scales((const double*)inits, scales);
	for (int i = 0; i < m->num_params(); i++) MSG("%.2f ", scales[i]); MSG("\n");
	//for (int i = 0; i < 5; i++) MSG("%.2f ", scales[i]); MSG("\n");
*/

/*
	// get samples
	HMC *sampler = new HMC(m, (const double *)inits, samples, deviance);
	sampler->sample(*Niter, *step_e, *step_L, *verbose);
*/

/*
	MH *sampler = new MH(m, (const double *)inits, samples, deviance);
	sampler->sample(*Niter, *thin, *step_e, true, *verbose);
*/

	SplineCovSampler *sampler = new SplineCovSampler(m);
	sampler->sample(samples, deviance, *Niter, (const double *)inits, *step_e, floor(*Niter/2), *verbose);

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
	double *lk, double *llik, double *lpri
) {

	// prior
	ModelSplineCov::Prior *model_prior = new ModelSplineCov::Prior();
	//model_prior->sd  = prior[0];
	//model_prior->var = pow(prior[0], 2);
	model_prior->mu_sd  = prior[0];
	model_prior->mu_var = pow(prior[0], 2.0);
	model_prior->tau_alpha = prior[1];
	model_prior->tau_beta  = prior[2];

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

	m->log_kernel((const double *)eval, lk, llik, lpri);

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
	//model_prior->sd  = prior[0];
	//model_prior->var = pow(prior[0], 2);
	model_prior->mu_sd  = prior[0];
	model_prior->mu_var = pow(prior[0], 2.0);
	model_prior->tau_alpha = prior[1];
	model_prior->tau_beta  = prior[2];

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

void spline_cov_scales(
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
	double *scales
) {

	// prior
	ModelSplineCov::Prior *model_prior = new ModelSplineCov::Prior();
	//model_prior->sd  = prior[0];
	//model_prior->var = pow(prior[0], 2);
	model_prior->mu_sd  = prior[0];
	model_prior->mu_var = pow(prior[0], 2.0);
	model_prior->tau_alpha = prior[1];
	model_prior->tau_beta  = prior[2];

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

	m->scales((const double*)eval, scales);

	delete m;
	delete model_prior;
	delete model_data;
}

} // end extern "C"
