#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "model.h"
#include "mh.h"
#include "utils.h"

#include <R.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>

// Specifies covariance mixture model
//
// Model of interest:
//     y_i   ind N(0, Sigma_i)
//
// Augmented model: (assume sigma^2 known and small)
//      y_i ind N(sum_l sqrt(w_il) alpha_il, sigma^2 I)
// alpha_il ind N(0, Omega_l)
//  Omega_l ind InvWish(nu, inv(S))

class ModelMixCov : public Model {
	friend class MixCovSampler;
public:
	class Prior {
	public:
		double        Omega_nu;
		const double *Omega_S;
	};

	class Data {
	public:
		int           n;             // number of observations
		int           k;             // number of sources
		const double *y;             // [n,k] observations
		int           L;             // number in mixture
		const int    *Nnz;           // [n] number of non-zero weights for this obs
		const int    *Mnz;           // [n,4] matrix of the non-zero indices
		const double *Wnz;           // [n,4] matrix of the non-zero weights
		double        sigma2;        // assumed known variance (should be small)
	};

	ModelMixCov(const Prior *prior, const Data *data);
	~ModelMixCov();

	//virtual bool log_kernel(const double *theta, double *lp, double *llik=NULL, double *lpri=NULL, bool do_lik=true, bool do_pri=true, int param=-1);

protected:
	const Prior *mPrior;
	const Data  *mData;
};

ModelMixCov::ModelMixCov(const Prior *prior, const Data *data) {
	mPrior = prior;
	mData  = data;
}

ModelMixCov::~ModelMixCov() {
}

int ModelMixCov::num_params() const {
	return(mNparams);
}

#if 0
bool ModelMixCov::log_kernel(const double *theta, double *lk, double *llik, double *lpri, bool do_lik, bool do_pri, int param) {
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
#endif

// do the sampling for this model
class MixCovSampler {
public:
	MixCovSampler(ModelMixCov *m);
	~MixCovSampler();

	bool sample(double *samples, double *deviance, int niter, const double *inits, double e, int nburn=0, bool verbose=false);

private:
	ModelMixCov *mModel;
	double *mAlpha;
	double *mOmega;
	int mIter;
};

MixCovSampler::MixCovSampler(ModelMixCov *m) {
	mModel   = m;
	mNparams = m->num_params();
	mNtheta  = m->num_theta();

	mPos = NULL;
	mPosCur = NULL;
}

MixCovSampler::~MixCovSampler() {
	free(mPos);
}

bool MixCovSampler::sample(double *samples, double *deviance,
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

void mix_cov_fit(
	// prior
	double *prior_nu, double *prior_S,
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
	ModelMixCov::Prior *model_prior = new ModelMixCov::Prior();
	model_prior->Omega_nu  = prior_nu[0];
	model_prior->Omega_S   = prior_S;

	// data
	ModelMixCov::Data *model_data = new ModelMixCov::Data();
	model_data->n   = n[0];
	model_data->k   = k[0];
	model_data->y   = (const double *)y;
	model_data->L   = L[0];
	model_data->Nnz = (const int *)Nnz;
	model_data->Mnz = (const int *)Mnz;
	model_data->Wnz = (const double *)Wnz;

	// model
	ModelMixCov *m = new ModelMixCov((const ModelMixCov::Prior *)model_prior, (const ModelMixCov::Data *)model_data);

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

	MixCovSampler *sampler = new MixCovSampler(m);
	sampler->sample(samples, deviance, *Niter, (const double *)inits, *step_e, floor(*Niter/2), *verbose);

	delete sampler;
	delete m;
	delete model_prior;
	delete model_data;
}

void mix_cov_lk(
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
	ModelMixCov::Prior *model_prior = new ModelMixCov::Prior();
	//model_prior->sd  = prior[0];
	//model_prior->var = pow(prior[0], 2);
	model_prior->mu_sd  = prior[0];
	model_prior->mu_var = pow(prior[0], 2.0);
	model_prior->tau_alpha = prior[1];
	model_prior->tau_beta  = prior[2];

	// data
	ModelMixCov::Data *model_data = new ModelMixCov::Data();
	model_data->n   = n[0];
	model_data->k   = k[0];
	model_data->y   = (const double *)y;
	model_data->L   = L[0];
	model_data->Nnz = (const int *)Nnz;
	model_data->Mnz = (const int *)Mnz;
	model_data->Wnz = (const double *)Wnz;

	// model
	Model *m = new ModelMixCov((const ModelMixCov::Prior *)model_prior, (const ModelMixCov::Data *)model_data);

	m->log_kernel((const double *)eval, lk, llik, lpri);

	delete m;
	delete model_prior;
	delete model_data;
}


} // end extern "C"
