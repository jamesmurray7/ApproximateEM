#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// -------- (gamma, eta) update --------
// Define the conditional expectation and then take Score AND Hessian 
// via forward differencing

// [[Rcpp::export]]
double Egammaeta(vec& gammaeta, mat& bmat, List S,
                 rowvec& K, mat& KK, mat& Fu, rowvec& Fi, vec& haz, int Delta, 
                 vec& w, vec& v){
  int gh = w.size();
  int nK = S.size();
  vec g = gammaeta.subvec(0, nK - 1);
  vec e = gammaeta.subvec(nK,nK + 1);
  // Rcout << "g: " << g << std::endl;
  // Rcout << "e: " << e << std::endl;
  vec tau = vec(Fu.n_rows);
  for(int i = 0; i < nK; i++){
    mat Si = S[i];
    tau += pow(g[i], 2.0) * diagvec(Fu * Si * Fu.t());
  }
  double rhs = 0.0;
  for(int l = 0; l < gh; l++){
    rhs += w[l] * as_scalar(haz.t() * exp(KK * e + Fu * (bmat.t() * g) + v[l] * pow(tau, 0.5)));
  }
  return as_scalar(Delta * (K * e + Fi * (bmat.t() * g)) - rhs);
}

// [[Rcpp::export]]
vec Sgammaeta(vec& gammaeta, mat& bmat, List S,
              rowvec& K, mat& KK, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v, double eps){
  int ge_size = gammaeta.size();
  vec out = vec(ge_size);
  double f0 = Egammaeta(gammaeta, bmat, S, K, KK, Fu, Fi, haz, Delta, w, v);
  for(int i = 0; i < ge_size; i++){
    vec ge = gammaeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammaeta[i] + xi * eps;
    double fdiff = Egammaeta(ge, bmat, S, K, KK, Fu, Fi, haz, Delta, w, v) - f0;
    out[i] = fdiff/(ge[i]-gammaeta[i]);
  }
  return out;
}

// [[Rcpp::export]]
mat Hgammaeta(vec& gammaeta, mat& bmat, List S,
              rowvec& K, mat& KK, mat& Fu, rowvec& Fi, vec& haz, int Delta, vec& w, vec& v, double eps){
  int ge_size = gammaeta.size();
  mat out = zeros<mat>(ge_size, ge_size);
  vec f0 = Sgammaeta(gammaeta, bmat, S, K, KK, Fu, Fi, haz, Delta, w, v, eps);
  for(int i = 0; i < ge_size; i++){
    vec ge = gammaeta;
    double xi = std::max(ge[i], 1.0);
    ge[i] = gammaeta[i] + xi * eps;
    vec fdiff = Sgammaeta(ge, bmat, S, K, KK, Fu, Fi, haz, Delta, w, v, eps) - f0;
    out.col(i) = fdiff/(ge[i]-gammaeta[i]);
  }
  return 0.5 * (out + out.t());
}
