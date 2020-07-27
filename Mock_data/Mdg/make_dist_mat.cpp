//sum.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix make_dist_mat(NumericMatrix X,
                     NumericMatrix X_obs){
  int n = X.rows();
  int n_obs = X_obs.rows();
  
  NumericMatrix dist_mat(n, n_obs);
  
  int n_covs = X.cols();
  
  double d=0;
  for(int i=0; i<n; i++){
    for(int j=0; j<n_obs; j++){
      d = 0;
      for(int k=0; k<n_covs; k++){
        d += pow(X(i, k) - X_obs(j, k), 2.0);
      }
      dist_mat(i, j) = d;
    }
  }
  return(dist_mat);
}