// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// The first to be used inside a loop (i.e. list indexing done in R)...
// [[Rcpp::export]]
double rcpp_e(const colvec& Yik, const mat& Xik, const mat& Zik, const vec& beta, const vec& b){
	colvec out = Yik - (Xik * beta + Zik * b);
	return as_scalar(out.t() * out);
}

// The next does the looping internally...
// OLD VERSION - NOT MV EXTENSION!!
/* NumericVector Ee(const Rcpp::List Y, const Rcpp::List X, const Rcpp::List Z, 
               const mat& beta, const Rcpp::List b, const int ids, const int K){
	mat M = zeros<mat>(ids, K);
	Rcpp::NumericVector e(K);
	for(int i = 0; i < ids; i++){
	  mat Yi = as<mat>(Y[i]);
	  mat bi = as<mat>(b[i]);
	  Rcpp::List Xi = X[i];
	  Rcpp::List Zi = Z[i];
	  
		for(int k = 0; k < K; k++){
		  colvec Yik = Yi.col(k);
		  mat Xik = as<mat>(Xi[k]);
		  mat Zik = as<mat>(Zi[k]);
		  colvec temp = Yik - (Xik * beta.row(k).t() + Zik * bi.row(k).t());
     	  M(i,k) = as_scalar(temp.t() * temp);
		}
	}
	// Work out column sums
	for(int j = 0; j < K; j++){
		e[j] = sum(M.col(j));
	}
	
	return e;
} */

// The next does the looping internally...
// [[Rcpp::export]]
NumericVector Ee(const Rcpp::List Y, const Rcpp::List X, const Rcpp::List Z, 
               const mat& beta, const Rcpp::List b, const Rcpp::List bbT, const int ids, const int K){
	mat M = zeros<mat>(ids, K);
	Rcpp::NumericVector e(K);
	for(int i = 0; i < ids; i++){
	  mat Yi = as<mat>(Y[i]);
	  mat bi = as<mat>(b[i]);
	  Rcpp::List Xi = X[i];
	  Rcpp::List Zi = Z[i];
	  Rcpp::List bbTi = bbT[i];
		for(int k = 0; k < K; k++){
		  colvec Yik = Yi.col(k);
		  mat Xik = as<mat>(Xi[k]);
		  mat Zik = as<mat>(Zi[k]);
		  mat bbTik = as<mat>(bbTi[k]);
		  colvec Resid = Yik - Xik * beta.row(k).t();
		  colvec temp = Yik - (Xik * beta.row(k).t() + Zik * bi.row(k).t());
     	  M(i,k) = as_scalar(
			Resid.t() * (Resid - 2.0 * (Zik * bi.row(k).t())) + trace(Zik.t() * Zik * bbTik)
		  );
		}
	}
	// Work out column sums
	for(int j = 0; j < K; j++){
		e[j] = sum(M.col(j));
	}
	
	return e;
}

// For E[b]
// [[Rcpp::export]]
List Eb(const List Y, const List X, const List Z, const List V, const mat& D, const vec & beta, const int ids){
	List b(ids);
	for(int i = 0; i < ids; i++){
		colvec Yi = Y[i];
		mat Xi = as<mat>(X[i]);
		mat Zi = as<mat>(Z[i]);
		mat Vi = as<mat>(V[i]);
		mat ZD = Zi * D;
		mat ZDZtV = Zi * D * Zi.t() + Vi;
		b[i] = ZD.t() * inv(ZDZtV) * (Yi - Xi * beta);
	}
	return b;
}

// VarCorr for b
// [[Rcpp::export]]
List covb(const List Z, const List V, const mat& Dinv, const int ids){
	List S(ids);
	for(int i = 0; i < ids; i++){
		mat Zi = as<mat>(Z[i]);
		mat Vi = as<mat>(V[i]);
		S[i] = inv(Zi.t() * inv(Vi) * Zi + Dinv);
	}
	return S;
}

// E[bbT]
// [[Rcpp::export]]
List EbbT(const List b, const List Sigma, const int ids){
	List bbT(ids);
	for(int i = 0; i < ids;  i++){
		colvec bi = b[i];
		mat Sigmai = as<mat>(Sigma[i]);
		bbT[i] = Sigmai + bi * bi.t();
	}
	return bbT;
}

// RHS for beta update
// [[Rcpp::export]]
List betaRHS(const List X, const List Y, const List Z, const List b, const int ids){
	List beta(ids);
	for(int i = 0; i < ids; i++){
		mat Xi = as<mat>(X[i]);
		mat Zi = as<mat>(Z[i]);
		colvec bi = b[i];
		colvec Yi = Y[i];
		beta[i] = Xi.t() * (Yi - Zi * bi);
	}
	return beta;
}
