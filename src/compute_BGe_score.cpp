#include <RcppArmadillo.h>
using namespace Rcpp;

#include<iostream>
#include<cmath>

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double compute_BGe_posterior(DataFrame df, 
                         unsigned alpha_w,
                         NumericVector nu,
                         double alpha_mu = 1.0
                         ) {
  
  double N = df.nrow();
  double n = df.ncol();
  double l = df.length();
  
  if (n == 0) return 0.0; // nothing to do on an empty variable set.
  
  // first term
  double res = l / 2.0 * (log(alpha_mu) - log(N + alpha_mu));
  
  // std::cout << "First step: " << res << std::endl;
  
  // StringVector nam = df.names();
  
  // Gamma_l ratio in the second term.
  arma::vec range = arma::linspace<arma::vec>(1, l, l);
  
  // std::cout << arma::lgamma((N + alpha_w - n + l + 1 - range) / 2) << std::endl;
  res += arma::sum(arma::lgamma((N + alpha_w - n + l + 1 - range) / 2) - 
    arma::lgamma((alpha_w - n + l + 1 - range) / 2));
  
  // std::cout << "Second step: " << res << std::endl;

  // leftover from the second term.
  res -= (l * N)/2 * log(M_PI);
  // std::cout << "Third step: " << res << std::endl;
          
  // third term, numerator.
  res += (alpha_w - n + l) / 2 * l * log((alpha_mu + alpha_w - n - 1)/(alpha_mu + 1));
  // std::cout << "Fourth step: " << res << std::endl;
  
  arma::mat m = as<arma::mat>(internal::convert_using_rfunction(df, "as.matrix"));
  
  arma::mat T(l, l, arma::fill::eye);
  T = T * (alpha_mu + alpha_w - n - 1)/(alpha_mu + 1);
  
  arma::vec v = arma::mean(m, 0).t() - as<arma::vec>(nu);
  arma::mat R = T + arma::cov(m) * (N - 1) + (N * alpha_w) / (N + alpha_w) * v * v.t();
  
  
  // std::cout << arma::cov(m, m) << std::endl;
  // std::cout << arma::cov(m) << std::endl;
  //std::cout << R << std::endl;
  
  res -= (N + alpha_w - n + l) / 2 * log(det(R));
  
  // third term, denominator. // Fix subsetting
//             nu = nu[vars]
//             local = data[, vars, drop = FALSE]
//               
//               res = res - (N + alpha_w - n + l)/2 * log(det(R))
//               
  return res;
}

// [[Rcpp::export]]
double compute_BGe_score(double N,
                         const arma::vec& means, 
                         const arma::mat& covmat,
                         List parents,
                         unsigned alpha_w,
                         const arma::vec& nu_vec,
                         double alpha_mu = 1.0
) {
  
  // double N = df.nrow();
  //double n = df.ncol();
  //double l = df.length();
  
  double n = means.size();
  
  // std::cout << "Number of variables: " << n << std::endl;
  
  if (n == 0) return 0.0; // nothing to do on an empty variable set.
  
  // StringVector df_names = df.names();
  
  // std::cout << df_names[0] << std::endl;
  
  double res = 0.0;
  
  for(int i = 0; i < n; ++i) {
    
    // std::cout << typeid(df_names[0]).name() << std::endl;
    // std::string node = as<std::string>(df_names[i]);
    // StringVector node_parents = parents[node];
    
    // arma::vec column = df[node];
    
    IntegerVector par = parents(i);
    double p = par.length();
    
    // root node
    if (p == 0) {
      
      // First term
      res += (log(alpha_mu) - log(N + alpha_mu)) / 2.0;
      // std::cout << "first term " << res << std::endl; 
      
      // Gamma_l ratio in the second term.
      res += lgamma((N + alpha_w - n + 1) / 2.0) -
        lgamma((alpha_w - n + 1) / 2.0);
      // std::cout << "second term 1 " << res << std::endl; 

      // leftover from the second term
      res -= N / 2.0 * log(M_PI);
      // std::cout << "second term 2 " << res << std::endl; 

      // third term, numerator
      
      // FIXED: Equations (19) and (20) in Geiger and Heckerman (2002)
      double t = alpha_mu * (alpha_w - n - 1) / (alpha_mu + 1);
      // double t = (alpha_mu + alpha_w - n - 1) / (alpha_mu + 1);
      
      res += (alpha_w - n + 1) / 2.0 * log(t);
      // std::cout << "Third term 1 " << res << std::endl; 
      
      // std::cout << typeid(arma::mean(column)).name() << std::endl;
      // std::cout << typeid(arma::var(column)).name() << std::endl;
      
      // double column_mean = arma::mean(column);
      // double column_var = arma::var(column);
      double nu = nu_vec(i);
        
      // third term, denominator
      // FIXED: Equation (17) in Geiger and Heckerman (2002)
      double r = t + covmat(i, i) * (N - 1) + (N * alpha_mu) / (N + alpha_mu) * (means(i) - nu) * (means(i) - nu);
      // double r = t + covmat(i, i) * (N - 1) + (N * alpha_w) / (N + alpha_w) * (means(i) - nu) * (means(i) - nu);
      
      res -= (N + alpha_w - n + 1) / 2.0 * log(r);
      
      // std::cout << "Third term 2 " << res << std::endl; 
      
    } else  { // p > 0
      
      // First term
      res += (log(alpha_mu) - log(N + alpha_mu)) / 2.0;
      
      // std::cout << "First term: " << res << std::endl;

      // Gamma_l ratio in the second term.
      res += lgamma((N + alpha_w - n + p + 1) / 2.0) - lgamma((alpha_w - n + p + 1) / 2.0);

      // leftover from the second term
      res -= N / 2.0 * log(M_PI);
      
      // std::cout << "Second term: " << res << std::endl;

      // third term, ratio of the determinants of the prior T matrices.
      // FIXED: Equations (19) and (20) in Geiger and Heckerman (2002)
      double t = alpha_mu * (alpha_w - n - 1) / (alpha_mu + 1);
      // double t = (alpha_mu + alpha_w - n - 1) / (alpha_mu + 1);
      res += (alpha_w - n + p + 1) / 2.0 * (p + 1) * log(t) -
              (alpha_w - n + p) / 2.0 * (p) * log(t);
      
      // std::cout << "Third term 1: " << res << std::endl;
    
      // std::vector<double> nu(p + 1);
      
      // std::vector<double> node_p(p + 1);
      
      
      // arma::uvec indices(p + 1);
      // indices << as<arma::uvec>(par);
      arma::uvec indices(p + 1);
      indices[0] = i;
      
      // TODO: simplify this
      for (int j = 0; j < p; ++j) indices[j + 1] = par[j];
      
// 
//       indices.resize(indices.size() + 1);
//       indices[indices.size() - 1] = i;
      
      // indices << i;
      
      // std::cout << indices;
      
      // indices << as<arma::uvec>(par);
      
      // node_p[0] = i;
      // 
      // for (int j = 0; j < p; ++j) {
      //   node_p[j+1] = par[j];
      // }
      // 
      // for (int j = 0; j <= p; ++j) {
      //   std::cout << node_p[j] << " ";
      // }
      // 
      // std::cout << std::endl;
      
      auto sm = covmat.submat(indices, indices);
      
      // std::cout << sm << std::endl;

      
      auto nu = nu_vec.elem(indices);
      // nu[0] = nu_vec[i];
      // 
      // for (int j = 0; j < p; ++j) {
      //   nu[j+1] = nu_vec[par[j]];
      // }
      
      //       
      // // third term, ratio of the determinants of the posterior R matrices.
      //       nu = nu[c(node, parents)]
      //       local = data[, c(node, parents), drop = FALSE]
      arma::vec v = means.elem(indices) - nu;
      
      // convoluted way to achieve T <- diag(t, p + 1, p + 1)
      arma::mat T(p+1, p+1, arma::fill::eye);
      T = T * t;
      
      // FIXED: Equation (17) in Geiger and Heckerman (2002)
      arma::mat R = T + sm * (N - 1) + (N * alpha_mu) / (N + alpha_mu) * v * v.t();
      // arma::mat R = T + sm * (N - 1) + (N * alpha_w) / (N + alpha_w) * v * v.t();
      
      
      arma::mat subR = R.submat(1, 1, p, p);
      
      // std::cout << R << std::endl;
      // std::cout << subR << std::endl;
      
      res += (N + alpha_w - n + p) / 2.0 * log(det(subR)) -
        (N + alpha_w - n + p + 1) / 2.0 * log(det(R));
      
      // std::cout << "Third term 2: " << res << std::endl;
        
        // diag(t, p + 1) + cov(local) * (N - 1) +
        //       (N * alpha_w) / (N + alpha_w) *
        //       outer(colMeans(local) - nu, colMeans(local) - nu)
        // 
        //       res = res + (N + alpha_w - n + p)/2 * log(det(R[-1, -1])) -
        //         (N + alpha_w - n + p + 1)/2 * log(det(R))
    }

  }
  
  // third term, denominator. // Fix subsetting
  //             nu = nu[vars]
  //             local = data[, vars, drop = FALSE]
  //               
  //               res = res - (N + alpha_w - n + l)/2 * log(det(R))
  //               
  return res;
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
  # library(bnlearn)
  # set.seed(1827)
  # dat <- data.frame(A = rnorm(100), B = rnorm(100), C = rnorm(100))
  # full_graph_1 = model2network("[A][B|A][C|A:B]")
  # parents_list <- lapply(full_graph_1$nodes, function(node) which(names(dat) %in% node$parents) - 1)
  # meanvec <- colMeans(dat)
  # covmat <- cov(dat)
  # compute_BGe_score(100, meanvec, covmat, parents_list, alpha_w = 5, nu = rep(0, 3))
  
  # compute_BGe_score(no_samples, ss$means, ss$covmat, parents_list, alpha_w = 5, nu_vec = rep(0, 3))
  
*/
