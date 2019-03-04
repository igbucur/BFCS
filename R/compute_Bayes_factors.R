#' Bayes Factor Computation for a three-variable submodel.
#'
#' @param corr Input 3x3 correlation matrix
#' @param no_samples Number of samples
#'
#' @return A vector containing the evidence of various models
#' @export
#'
#' @examples
compute_Bayes_factors <- function(corr, no_samples) {
  stopifnot(nrow(corr) == 3)

  # if (eigen(corr)$values[3] < 0) return (NA)

  nu <- nrow(corr) + 1 # uniform
  Nnu <- (no_samples + nu) / 2
  
  det_cor <- det(corr)

  # c_1(n, v)
  c1 <- (no_samples + nu - 2) / (nu - 2)
  
  # c_2(n, v)
  # We have to use logarithms for stability when Nnu is large
  c2 <- exp(lgamma(Nnu) - lgamma(Nnu - 0.5)) * gamma((nu - 1) / 2) / gamma(nu / 2)
  # c2 <- gamma(Nnu) * gamma((nu - 1) / 2) / (gamma(Nnu - 0.5) * gamma(nu / 2))
  
  # Approximation to c2 is not significantly faster, but more stable
  # c2 <- sqrt((no_samples + nu - 1.5) / (nu - 1.5))
  # c2 <- gamma((N + nu) / 2) * gamma((nu - 1) / 2) /
  #   (gamma((N + nu - 1) / 2) * gamma(nu / 2))

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


#' Vectorized compute Bayes factors of covariance structures.
#'
#' @param c12 vector of correlations between X_1 and X_2
#' @param c13 vector of correlations between X_1 and X_3
#' @param c23 vector of correlations between X_2 and X_3
#' @param no_samples number of samples
#'
#' @return Bayes factors for all correlations in the vectors
#' @export
#'
#' @examples
compute_Bayes_factors_vectorized <- function(c12, c13, c23, no_samples) {
  
  vecl <- length(c12)
  stopifnot(vecl == length(c13))
  stopifnot(vecl == length(c23))
  
  # if (eigen(corr)$values[3] < 0) return (NA)
  
  nu <- 4 # uniform
  Nnu <- (no_samples + nu) / 2
  
  # c_1(n, v)
  c1 <- (no_samples + nu - 2) / (nu - 2)
  
  # c_2(n, v)
  # We have to use logarithms for stability when Nnu is large
  c2 <- exp(lgamma(Nnu) - lgamma(Nnu - 0.5)) * gamma((nu - 1) / 2) / gamma(nu / 2)
  # c2 <- gamma(Nnu) * gamma((nu - 1) / 2) / (gamma(Nnu - 0.5) * gamma(nu / 2))
  
  # Approximation to c2 is not significantly faster, but more stable
  # c2 <- sqrt((no_samples + nu - 1.5) / (nu - 1.5))
  # c2 <- gamma((N + nu) / 2) * gamma((nu - 1) / 2) /
  #   (gamma((N + nu - 1) / 2) * gamma(nu / 2))
  
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



yeast_BFCM_prob <- function(L_i, T_i, T_j, prior_structure = NULL) {
  data(yeast)
  
  if (is.null(prior_structure)) {
    prior_structure <- uniform_prior_GRN_triplet()
  }
  
  Bf <- compute_Bayes_factors(cor(cbind(L_i, T_i, T_j)), 112) * prior_structure
  Bf / sum(Bf)
}


