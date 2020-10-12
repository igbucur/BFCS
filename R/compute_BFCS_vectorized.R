#' Function for computing the Bayes Factors of Covariance Structures for a 
#' three-variable causal model.
#'
#' @param corr Numeric 3x3 correlation matrix.
#' @param num_samples Integer number of samples.
#'
#' @return A vector containing the Bayes factor of all the eleven possible
#' conditional independence models on three variables relative to the fully
#' connected model.
#' @export
#'
#' @examples
#' compute_Bayes_factors(diag(3), 100)
#' 
compute_Bayes_factors <- function(corr, num_samples) {
  
  if (nrow(corr) != 3 || ncol(corr) != 3) {
    stop("Input correlation matrix has wrong dimensions, should be 3x3.")
  }

  nu <- nrow(corr) + 1 # corresponds to uniform off-diagonal correlation
  Nnu <- (num_samples + nu) / 2
  
  det_cor <- det(corr)

  # Compute c_1(n, v)
  c1 <- (num_samples + nu - 2) / (nu - 2)
  
  # Compute c_2(n, v) - We have to use logarithms for stability when Nnu is large
  c2 <- exp(lgamma(Nnu) - lgamma(Nnu - 0.5)) * gamma((nu - 1) / 2) / gamma(nu / 2)

  f <- function(i, j) (1 - corr[i, j]^2)

  ## 1. 1 \indep 2 \indep 3 -- 1 DAG
  BF_all_indep <- c1 * c2 * det_cor^Nnu

  ## 2. INDEPENDENT -- 6 DAG
  
  # 1 \indep (2, 3) -- 2 DAG
  BF_indep_1 <- c1 * (det_cor / f(2, 3))^Nnu

  # 2 \indep (3, 1) -- 2 DAG
  BF_indep_2 <- c1 * (det_cor / f(3, 1))^Nnu

  # 3 \indep (1, 2) -- 2 DAG
  BF_indep_3 <- c1 * (det_cor / f(1, 2))^Nnu

  ## 3. CAUSAL -- 9 DAG

  # 2 \indep 3 \given 1 -- 3 DAG
  BF_causal_123 <- c2 * (det_cor / (f(1, 2) * f(1, 3)))^Nnu
  
  # 3 \indep 1 \given 2 -- 3 DAG
  BF_causal_231 <- c2 * (det_cor / (f(2, 3) * f(2, 1)))^Nnu

  # 1 \indep 2 \given 3 -- 3 DAG
  BF_causal_312 <- c2 * (det_cor / (f(3, 1) * f(3, 2)))^Nnu

  ## 4. ACAUSAL -- 3 DAG

  # 1 \indep 2 -- 1 DAG
  BF_acausal_12 <- c1 * (f(1, 2)^(Nnu - 0.5)) / c2

  # 2 \indep 3 -- 1 DAG
  BF_acausal_23 <- c1 * (f(2, 3)^(Nnu - 0.5)) / c2
  
  # 3 \indep 1 -- 1 DAG
  BF_acausal_31 <- c1 * (f(3, 1)^(Nnu - 0.5)) / c2

  ## 5. FULL -- 6 DAG
  BF_no_indep <- 1 

  factors <- c(
    BF_all_indep,
    BF_indep_1, BF_indep_2, BF_indep_3,
    BF_causal_123, BF_causal_231, BF_causal_312,
    BF_acausal_12, BF_acausal_23, BF_acausal_31,
    BF_no_indep
  )

  factors
}


#' Vectorized function for efficiently computing the Bayes Factors of many 
#' three-variable covariance structures at the same time.
#'
#' @param c12 Vector of correlations between X_1 and X_2.
#' @param c13 Vector of correlations between X_1 and X_3.
#' @param c23 Vector of correlations between X_2 and X_3.
#' @param num_samples Integer number of samples.
#'
#' @return BFCS for all correlation combinations contained in the vectors
#' @export
#'
#' @examples
#' cor_matrices <- apply(rWishart(100, 4, diag(3)), 3, cov2cor)
#' compute_Bayes_factors_vectorized(
#'   c12 = cor_matrices[2, ], 
#'   c13 = cor_matrices[3, ], 
#'   c23 = cor_matrices[6, ], 
#'   num_samples = 1000
#'  )
compute_Bayes_factors_vectorized <- function(c12, c13, c23, num_samples) {
  
  vecl <- length(c12)
  stopifnot(vecl == length(c13))
  stopifnot(vecl == length(c23))
  
  nu <- 4 # uniform
  Nnu <- (num_samples + nu) / 2
  
  # c_1(n, v)
  c1 <- (num_samples + nu - 2) / (nu - 2)
  
  # Compute c_2(n, v) - We have to use logarithms for stability when Nnu is large
  c2 <- exp(lgamma(Nnu) - lgamma(Nnu - 0.5)) * gamma((nu - 1) / 2) / gamma(nu / 2)
  
  f23 <- 1 - c23 * c23
  f13 <- 1 - c13 * c13
  f12 <- 1 - c12 * c12
  
  det_cor <- 1 + 2 * c12 * c13 * c23 - c12 * c12 - c23 * c23 - c13 * c13
  
  Bf <- matrix(0, vecl, 11)
  colnames(Bf) <- c("empty", "indep_1", "indep_2", "indep_3", 
                 "causal_123", "causal_231", "causal_312",
                 "acausal_12", "acausal_23", "acausal_31", "full")
  
  
  ## 1. 1 \indep 2 \indep 3 -- 1 DAG
  Bf[, 1] <- c1 * c2 * det_cor^Nnu
  
  ## 2. INDEPENDENT -- 6 DAG
  
  # 1 \indep (2, 3) -- 2 DAG
  Bf[, 2] <- c1 * (det_cor / f23)^Nnu
  
  # 2 \indep (3, 1) -- 2 DAG
  Bf[, 3] <- c1 * (det_cor / f13)^Nnu
  
  # 3 \indep (1, 2) -- 2 DAG
  Bf[, 4] <- c1 * (det_cor / f12)^Nnu
  
  ## 3. CAUSAL -- 9 DAG
  
  # 2 \indep 3 \given 1 -- 3 DAG
  Bf[, 5] <- c2 * (det_cor / (f12 * f13))^Nnu
  
  # 3 \indep 1 \given 2 -- 3 DAG
  Bf[, 6] <- c2 * (det_cor / (f23 * f12))^Nnu
  
  # 1 \indep 2 \given 3 -- 3 DAG
  Bf[, 7] <- c2 * (det_cor / (f13 * f23))^Nnu
  
  ## 4. ACAUSAL -- 3 DAG
  
  # 1 \indep 2 -- 1 DAG
  Bf[, 8] <- c1 * (f12^(Nnu - 0.5)) / c2
  
  # 2 \indep 3 -- 1 DAG
  Bf[, 9] <- c1 * (f23^(Nnu - 0.5)) / c2
  
  # 3 \indep 1 -- 1 DAG
  Bf[, 10] <- c1 * (f13^(Nnu - 0.5)) / c2
  
  ## 5. FULL -- 6 DAG
  Bf[, 11] <- 1 
  
  Bf
}

