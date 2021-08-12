#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Split out getxi
// [[Rcpp::export]]
colvec getxi(const rowvec& tausurv, const colvec& musurv, double v, const rowvec& haz){
  rowvec xi = haz % (musurv.t() % exp(v * tausurv));
  return xi.t();
}

// [[Rcpp::export]]
mat gammaCalc(const vec& gamma,const mat& tautilde, const rowvec& tausurv, const rowvec& tau2surv, const colvec& musurv, 
              const vec& w, const vec& v, const mat& Fu, const rowvec& haz, const mat& bb, int L, int gh){
	mat M = zeros<mat>(L, L);
	for(int i = 0; i < L; i++){
		rowvec bM = bb.row(i);
		for(int j = 0; j < L; j++){
			rowvec bL = bb.row(j);
			// Rcpp::Rcout << " (i,j)= " << i << "," << j << std::endl; 
			if(i!=j){ // The cross-terms...
				for(int k = 0; k < gh; k++){
					const colvec xi = getxi(tausurv, musurv, v(k), haz);
					// Rcpp::Rcout << "xi: " << xi << std::endl;
					const colvec temp1 = xi % tausurv.t() % xi;
					// Rcpp::Rcout << "xi * tau.surv  * xi: " << temp1 << std::endl;
					const colvec temp2 = temp1 % tau2surv.t();
					// Rcpp::Rcout << "xi * tau.surv * xi * tau2surv: " << temp2 << std::endl;
					const colvec temp3 = xi % xi % tau2surv.t();
					// Rcpp::Rcout << "xi * xi * tau2surv: " << temp3 << std::endl;
					
					// Rcpp::Rcout << "Whole first line: " << w(k) * bL* Fu.t() * (xi % (Fu * bM.t())) << std::endl;
					// Rcpp::Rcout << "Whole second line: " << gamma[i] * w[k] * v[k] * bL * Fu.t() * (xi % tau2surv.t() % tautilde.row(i).t()) << std::endl;
					// Rcpp::Rcout << "Whole third line: " << 2 * gamma[j] * w[k] * v[k] * temp1.t() * (Fu * bM.t()) << std::endl;
					// Rcpp::Rcout << "Whole fourth line: " << 2 * gamma[i] * gamma[j] * w[k] * v[k] * v[k] * temp2.t() * tautilde.row(i).t() << std::endl;
					// Rcpp::Rcout << "Whole fifth line: " << gamma[i] * gamma[j] * v[k] * w[k] * temp3.t() * tautilde.row(i).t() << std::endl;
					
					M(i,j) += as_scalar(w(k) * bL * Fu.t() * (xi % (Fu * bM.t())) + 
					          gamma[i] * w[k] * v[k] * bL * Fu.t() * (xi % tau2surv.t() % tautilde.row(i).t()) +
					          2 * gamma[j] * w[k] * v[k] * temp1.t() * (Fu * bM.t()) + 
					          2 * gamma[i] * gamma[j] * w[k] * v[k] * v[k] * temp2.t() * tautilde.row(i).t() + 
					          gamma[i] * gamma[j] * v[k] * w[k] * temp3.t() * tautilde.row(i).t());
				}
			}else{ // The diagonal terms d2f/dgdg...
				for(int k = 0; k < gh; k++){
					const colvec xi = getxi(tausurv, musurv, v(k), haz);
					const colvec xitau = xi % tausurv.t();
					const colvec temp1 = xi % tausurv.t() % xi;
					const colvec temp2 = temp1 % tau2surv.t();
					const colvec temp3 = xi % xi % tau2surv.t();
					
					M(i,j) += as_scalar(w(k) * bL * Fu.t() * (xi % (Fu * bM.t())) + 
						      gamma[i] * w[k] * v[k] * bL * Fu.t() * (xi % tau2surv.t() % tautilde.row(i).t()) + 
						      v[k] * w[k] * xitau.t() * xi + 
						      2 * gamma[i] * v[k] * w[k] * temp1.t() * (Fu * bL.t()) + 
						      2 * gamma[i] * gamma[i] * w[k] * v[k] * v[k] * temp2.t() * tautilde.row(i).t() + 
						      gamma[i] * gamma[i] * v[k] * w[k] * temp3.t() * tautilde.row(i).t());
					
				}
			}
		}
	}
	return M;
}

// [[Rcpp::export]]
mat gamma2Calc(const vec& gamma,const mat& tautilde, const rowvec& tausurv, const rowvec& tau2surv, const colvec& musurv, 
               const vec& w, const vec& v, const mat& Fu, const rowvec& haz, const mat& bb, int L, int gh){
	mat M = zeros<mat>(L, L);
	for(int i = 0; i < L; i++){
		rowvec bM = bb.row(i);
		for(int j = 0; j < L; j++){
			rowvec bL = bb.row(j);
			for(int k = 0; k < gh; k++){
				const colvec xi = getxi(tausurv, musurv, v(k), haz);
				// Rcpp::Rcout << "xi: " << xi << std::endl;
				const colvec temp1 = xi % tausurv.t() % xi;
				// Rcpp::Rcout << "xi * tau.surv  * xi: " << temp1 << std::endl;
				const colvec temp2 = temp1 % tau2surv.t();
				// Rcpp::Rcout << "xi * tau.surv * xi * tau2surv: " << temp2 << std::endl;
				const colvec temp3 = xi % xi % tau2surv.t();
				// Rcpp::Rcout << "xi * xi * tau2surv: " << temp3 << std::endl;
				const colvec xitau = xi % tausurv.t();
				
				// Rcpp::Rcout << "Whole first line: " << w(k) * bL* Fu.t() * (xi % (Fu * bM.t())) << std::endl;
				// Rcpp::Rcout << "Whole second line: " << gamma[i] * w[k] * v[k] * bL * Fu.t() * (xi % tau2surv.t() % tautilde.row(i).t()) << std::endl;
				// Rcpp::Rcout << "Whole third line: " << 2 * gamma[j] * w[k] * v[k] * temp1.t() * (Fu * bM.t()) << std::endl;
				// Rcpp::Rcout << "Whole fourth line: " << 2 * gamma[i] * gamma[j] * w[k] * v[k] * v[k] * temp2.t() * tautilde.row(i).t() << std::endl;
				// Rcpp::Rcout << "Whole fifth line: " << gamma[i] * gamma[j] * v[k] * w[k] * temp3.t() * tautilde.row(i).t() << std::endl;
				if(i<j){ // The upper triangle
					M(i,j) += as_scalar(w(k) * bL * Fu.t() * (xi % (Fu * bM.t())) + 
						  gamma[i] * w[k] * v[k] * bL * Fu.t() * (xi % tau2surv.t() % tautilde.row(i).t()) +
						  2 * gamma[j] * w[k] * v[k] * temp1.t() * (Fu * bM.t()) + 
						  2 * gamma[i] * gamma[j] * w[k] * v[k] * v[k] * temp2.t() * tautilde.row(i).t() + 
						  gamma[i] * gamma[j] * v[k] * w[k] * temp3.t() * tautilde.row(i).t());
				}else if(i == j){ // The diagonal terms
					M(i,j) += as_scalar(w(k) * bL * Fu.t() * (xi % (Fu * bM.t())) + 
						      gamma[i] * w[k] * v[k] * bL * Fu.t() * (xi % tau2surv.t() % tautilde.row(i).t()) + 
						      v[k] * w[k] * xitau.t() * xi + 
						      2 * gamma[i] * v[k] * w[k] * temp1.t() * (Fu * bL.t()) + 
						      2 * gamma[i] * gamma[i] * w[k] * v[k] * v[k] * temp2.t() * tautilde.row(i).t() + 
						      gamma[i] * gamma[i] * v[k] * w[k] * temp3.t() * tautilde.row(i).t());
				}
			}
		M(j,i) = M(i,j);
		}
	}
	return M;
}
