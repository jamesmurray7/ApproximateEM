#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double bllC(vec& b, const colvec& Y, const mat& X, const mat& Z, const mat& V, const mat& D,
		   int mi, const rowvec& K, const int Delta, const double l0i, const rowvec& Fi,
		   const rowvec& l0u, const mat& Fu, const vec& g, const vec& beta, const vec& eta,
		   const vec& gr, const rowvec& rvFi){
	double temp = 0.0;
	if(Delta == 1) temp = log(l0i);
	colvec resid = Y - X * beta - Z * b;
	return -1.0 * as_scalar(-mi/2.0 * log(2.0 * M_PI) - 0.5 * log(det(V)) -0.5 * resid.t() * V.i() * resid +
	                 -5.0 * log(2.0 * M_PI) - 0.5 * log(det(D)) - 0.5 * b.t() * D.i() * b + 
	                 temp + Delta * (K * eta + rvFi * (gr % b)) - l0u * (kron(exp(K * eta), exp(repmat(Fu, 1, 5) * (gr % b)))));
}

// Gradient function
// [[Rcpp::export]]
colvec gradC(vec& b, const colvec& Y, const mat& X, const mat& Z, const mat& V, const mat& D,
		     int mi, const rowvec& K, const int Delta, const double l0i, const rowvec& Fi,
		     const rowvec& l0u, const mat& Fu, const vec& g, const vec& beta, const vec& eta,
		     const vec& gr, const rowvec& rvFi){
	colvec resid = Y - X * beta - Z * b;
	return -1.0 * (Z.t() * V.i() * resid - D.i() * b + Delta * rvFi.t() % gr + 
	                 - repmat(Fu, 1, 5).t() * (l0u.t() % (kron(exp(K * eta), exp(repmat(Fu, 1, 5) * (gr % b))))) % gr);

}
