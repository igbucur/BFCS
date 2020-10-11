#include <RcppArmadillo.h>
using namespace Rcpp;

#include<iostream>
#include<cmath>


// [[Rcpp::depends(RcppArmadillo)]]


//' Compute Bayesian Gaussian equivalent (BGe) score for a given Bayesian network
//' on given sufficient statistics from data.
//' 
//' @details For a better understanding of how the BGe score is obtained see the article 
//' [Geiger and Heckerman (2002)](https://projecteuclid.org/euclid.aos/1035844981)
//' as well as the discussion in [Kuipers et al. (2014)](https://projecteuclid.org/euclid.aos/1407420013).
//'
//' @param N Integer number of observations.
//' @param means Numeric vector containing variable means.
//' @param covmat Numeric covariance matrix.
//' @param parents List of parents for each node to describe Bayesian network.
//' @param alpha_w Numeric parameter of the covariance's Wishart prior (scale matrix).
//' @param nu_vec Numeric parameter of the covariance's Wishart prior (degrees of freedom)
//' @param alpha_mu Numeric parameter of the mean's Gaussian prior (expected value).
//'
//' @return Bayesian Gaussian equivalent (BGe) score for given data on the
//' given Bayesian network.
// [[Rcpp::export]]
double compute_BGe_score(double N,
                         const arma::vec& means, 
                         const arma::mat& covmat,
                         List parents,
                         unsigned alpha_w,
                         const arma::vec& nu_vec,
                         double alpha_mu = 1.0
) {
  
  // Get number of variables
  double n = means.size();
  if (n == 0) return 0.0; // nothing to do on an empty variable set.
  
  // Initialize score variable
  double BGe_score = 0.0;
  
  // For each variable
  for(int i = 0; i < n; ++i) {
    
    // Get parents of current variable and their number
    IntegerVector par = parents(i);
    double p = par.length();
    
    // no parents means a root node
    if (p == 0) {
      
      // Add first term from BGe expression
      BGe_score += (log(alpha_mu) - log(N + alpha_mu)) / 2.0;
      
      // Add second term from BGe expression, Gamma_l ratio
      BGe_score += lgamma((N + alpha_w - n + 1) / 2.0) -
        lgamma((alpha_w - n + 1) / 2.0);
      // leftover from the second term
      BGe_score -= N / 2.0 * log(M_PI);
      
      // third term, numerator
      
      // FIXED: Equations (19) and (20) in Geiger and Heckerman (2002)
      double t = alpha_mu * (alpha_w - n - 1) / (alpha_mu + 1);
      
      BGe_score += (alpha_w - n + 1) / 2.0 * log(t);

      // third term, denominator
      double nu = nu_vec(i);
        
      // FIXED: Equation (17) in Geiger and Heckerman (2002)
      double r = t + covmat(i, i) * (N - 1) + (N * alpha_mu) / (N + alpha_mu) * (means(i) - nu) * (means(i) - nu);

      BGe_score -= (N + alpha_w - n + 1) / 2.0 * log(r);
      
    } else  { // p > 0
      
      // Add first term from BGe expression
      BGe_score += (log(alpha_mu) - log(N + alpha_mu)) / 2.0;

      // Add second term from BGe expression, Gamma_l ratio
      BGe_score += lgamma((N + alpha_w - n + p + 1) / 2.0) - lgamma((alpha_w - n + p + 1) / 2.0);
      // leftover from the second term
      BGe_score -= N / 2.0 * log(M_PI);
      
      // Add third term, ratio of the determinants of the prior T matrices.
      // FIXED: Equations (19) and (20) in Geiger and Heckerman (2002)
      double t = alpha_mu * (alpha_w - n - 1) / (alpha_mu + 1);
      
      BGe_score += (alpha_w - n + p + 1) / 2.0 * (p + 1) * log(t) -
              (alpha_w - n + p) / 2.0 * (p) * log(t);
      
      // Add third term, ratio of the determinants of the posterior R matrices.
      
      arma::uvec indices(p + 1); indices[0] = i;
      for (int j = 0; j < p; ++j) indices[j + 1] = par[j];
      
      auto sm = covmat.submat(indices, indices); // extract submatrix of covmat
      
      auto nu = nu_vec.elem(indices);
      arma::vec v = means.elem(indices) - nu;
      
      // convoluted way of achieving T <- diag(t, p + 1, p + 1)
      arma::mat T(p+1, p+1, arma::fill::eye); T = T * t;
      
      // FIXED: Equation (17) in Geiger and Heckerman (2002)
      arma::mat R = T + sm * (N - 1) + (N * alpha_mu) / (N + alpha_mu) * v * v.t();
      arma::mat subR = R.submat(1, 1, p, p); // extract submatrix of R
      
      BGe_score += (N + alpha_w - n + p) / 2.0 * log(det(subR)) -
        (N + alpha_w - n + p + 1) / 2.0 * log(det(R));

    }

  }
  
  return BGe_score;
}

