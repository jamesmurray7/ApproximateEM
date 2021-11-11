#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// S is this -->
// lapply(Sigmai.store, function(x) lapply(split(seq(6), rep(1:3, each = 2)), function(y) x[y,y]))

// Lambda update, does all AFTER the E-step 
// [[Rcpp::export]]
mat lambdaUpdate(const List survtimes, const vec& ft,
				 const vec& gamma, const vec& eta, 
				 const List K, const List S, List b,
				 const int id, 
				 const vec& w, const vec& v, const int nodes, const int nK){
	mat store = zeros<mat>(ft.size(), id); // Initialise the matrix
	// Start loop over i subjects
	for(int i = 0; i < id; i++){
		vec survtimes_i = survtimes[i];    // This id's survived time indices   
		List Si = S[i];
		List bi = b[i];
		rowvec Ki = K[i];           // Start loop over subject i's j survived times     
		for(int j = 0; j < survtimes_i.size(); j++){
			rowvec Fst = NumericVector::create(1.0, ft[j]);
			double tau = 0.0;
			vec rhs = NumericVector::create(0.0, 0.0); // intslope hardcode...
			
			for(int k = 0; k < nK; k++){   // Loop over the nK longitudinal responses
				mat Sk = as<mat>(Si[k]);
				vec bk = bi[k];
				tau += as_scalar(pow(gamma[k], 2.0) * Fst * Sk * Fst.t());
				rhs += gamma[k] * bk;
			}
			double mu = as_scalar(exp(Ki * eta + Fst * rhs));
			// Rcpp::Rcout << "mu = " << mu << std::endl;								  
			for(int l = 0; l < nodes; l++){ // Finally loop over gh nodes
				store(j,i) += as_scalar(w[l] * mu * exp(v[l] * sqrt(tau)));
			}
		}
	}
	
	return store;
}
