#include <TMB.hpp>
#include <math.h>
#include <stdio.h>

const double pi = 3.141592653589793238462643383279502884;


template<class Type>
  Type objective_function<Type>::operator() ()
{
  
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  DATA_MATRIX(X);
  DATA_VECTOR(Y);
  DATA_VECTOR(pops);
  DATA_SPARSE_MATRIX(A);
  DATA_STRUCT(spde,spde_t);
  DATA_SCALAR(prior_rho_min);
  DATA_SCALAR(prior_rho_prob);
  DATA_SCALAR(prior_sigma_max);
  DATA_SCALAR(prior_sigma_prob);
  DATA_SCALAR(nu);
  
  PARAMETER(beta_0);
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(S);
  PARAMETER(log_rho);
  PARAMETER(log_sigma);
  
  Type f=0;
  Type sigma = exp(log_sigma);
  Type rho = exp(log_rho);
  Type kappa = sqrt(8.0*nu) / rho;
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  
  Type lambdatilde1 = -log(prior_rho_prob) * prior_rho_min;
  Type lambdatilde2 = -log(prior_sigma_prob) / prior_sigma_max;
  Type pcdensity = lambdatilde1 * lambdatilde2 * pow(rho, -2) * exp(-lambdatilde1 * pow(rho, -1) - lambdatilde2 * sigma);
  f -= log(pcdensity) + log_rho + log_sigma;
  Type scaling_factor = sqrt(exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * pi * pow(kappa, 2*nu)));

  //priors
  Type beta_mean = 0.0;
  Type beta_sd = 1.0;
  int n_beta = beta.size();
  int n_obs = Y.size();
  
  for(int i=0; i<n_beta; i++){
    f -= dnorm(beta[i], beta_mean, beta_sd, true);
  }
  
  vector<Type> field = A * S;
  f += SCALE(GMRF(Q), sigma / scaling_factor)(S);

  vector<Type> log_rate = beta_0 + (X * beta) + field;

  for(int i=0; i<n_obs; i++){
    f -= dbinom_robust(Y(i), pops(i), log_rate(i), true);
  }
  
  return(f);
  }
