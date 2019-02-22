#' Uniform prior on three-variable DAG structures
#'
#' @return Logical vector indicating the probability of ...
#' @export
#'
#' @examples
uniform_prior_DAG_3var <- function() {
  
  n_all_indep <- 1 # 1 \indep 2 \indep 3
  
  n_indep_1 <- 2 # 1 \indep (2, 3)
  n_indep_2 <- 2 # 2 \indep (3, 1)
  n_indep_3 <- 2 # 3 \indep (1, 2)
  
  n_causal_123 <- 3 # 2 \indep 3 \given 1
  n_causal_231 <- 3 # 3 \indep 1 \given 2
  n_causal_312 <- 3 # 1 \indep 2 \given 3

  n_acausal_12 <- 1 # 1 \indep 2
  n_acausal_23 <- 1 # 2 \indep 3
  n_acausal_31 <- 1 # 3 \indep 1
  
  n_no_indep <- 6 # all dependent
  
  structures <- c(
    n_all_indep,
    n_indep_1, n_indep_2, n_indep_3,
    n_causal_123, n_causal_231, n_causal_312,
    n_acausal_12, n_acausal_23, n_acausal_31,
    n_no_indep
  )
  
  structures / sum(structures)
}

#' Uniform prior on DAG structures (for GRN)
#'
#' @return
#' @export
#'
#' @examples
uniform_prior_GRN_DAG <- function() {
  
  # 1 == L_i, 2 == T_i, 3 == T_j
  
  n_all_indep <- 1 # 1 \indep 2 \indep 3
  
  n_indep_1 <- 2 # 1 \indep (2, 3)
  n_indep_2 <- 1 # 2 \indep (3, 1)
  n_indep_3 <- 1 # 3 \indep (1, 2)
  
  n_causal_123 <- 1 # 2 \indep 3 \given 1
  n_causal_231 <- 1 # 3 \indep 1 \given 2
  n_causal_312 <- 1 # 1 \indep 2 \given 3
  
  n_acausal_12 <- 1 # 1 \indep 2
  n_acausal_23 <- 0 # 2 \indep 3
  n_acausal_31 <- 1 # 3 \indep 1
  
  n_no_indep <- 2 # all dependent
  
  structures <- c(
    n_all_indep,
    n_indep_1, n_indep_2, n_indep_3,
    n_causal_123, n_causal_231, n_causal_312,
    n_acausal_12, n_acausal_23, n_acausal_31,
    n_no_indep
  )
  
  structures / sum(structures)
}

#' Uniform prior on DMAG structures
#'
#' @return
#' @export
#'
#' @examples
uniform_prior_GRN_DMAG <- function() {
  
  # 1 == L_i, 2 == T_i, 3 == T_j
  
  n_all_indep <- 1 # 1 \indep 2 \indep 3
  
  n_indep_1 <- 3 # 1 \indep (2, 3)
  n_indep_2 <- 1 # 2 \indep (3, 1)
  n_indep_3 <- 1 # 3 \indep (1, 2)
  
  n_causal_123 <- 1 # 2 \indep 3 \given 1
  n_causal_231 <- 1 # 3 \indep 1 \given 2
  n_causal_312 <- 1 # 1 \indep 2 \given 3
  
  n_acausal_12 <- 2 # 1 \indep 2
  n_acausal_23 <- 0 # 2 \indep 3
  n_acausal_31 <- 2 # 3 \indep 1
  
  n_no_indep <- 3 # all dependent
  
  structures <- c(
    n_all_indep,
    n_indep_1, n_indep_2, n_indep_3,
    n_causal_123, n_causal_231, n_causal_312,
    n_acausal_12, n_acausal_23, n_acausal_31,
    n_no_indep
  )
  
  structures / sum(structures)
}

#' Uniform prior on MAG structures
#'
#' @return
#' @export
#'
#' @examples
uniform_prior_GRN_MAG <- function() {
  
  # 1 == L_i, 2 == T_i, 3 == T_j
  
  n_all_indep <- 1 # 1 \indep 2 \indep 3
  
  n_indep_1 <- 4 # 1 \indep (2, 3)
  n_indep_2 <- 2 # 2 \indep (3, 1)
  n_indep_3 <- 2 # 3 \indep (1, 2)
  
  n_causal_123 <- 4 # 2 \indep 3 \given 1
  n_causal_231 <- 3 # 3 \indep 1 \given 2
  n_causal_312 <- 3 # 1 \indep 2 \given 3
  
  n_acausal_12 <- 2 # 1 \indep 2
  n_acausal_23 <- 0 # 2 \indep 3
  n_acausal_31 <- 2 # 3 \indep 1
  
  n_no_indep <- 6 # all dependent
  
  structures <- c(
    n_all_indep,
    n_indep_1, n_indep_2, n_indep_3,
    n_causal_123, n_causal_231, n_causal_312,
    n_acausal_12, n_acausal_23, n_acausal_31,
    n_no_indep
  )
  
  structures / sum(structures)
}


#' Uniform prior on ADMG structures
#'
#' @return
#' @export
#'
#' @examples
uniform_prior_GRN_ADMG <- function() {
  
  # 1 == L_i, 2 == T_i, 3 == T_j
  
  n_all_indep <- 1 # 1 \indep 2 \indep 3
  
  n_indep_1 <- 5 # 1 \indep (2, 3)
  n_indep_2 <- 1 # 2 \indep (3, 1)
  n_indep_3 <- 1 # 3 \indep (1, 2)
  
  n_causal_123 <- 1 # 2 \indep 3 \given 1
  n_causal_231 <- 1 # 3 \indep 1 \given 2
  n_causal_312 <- 1 # 1 \indep 2 \given 3
  
  n_acausal_12 <- 3 # 1 \indep 2
  n_acausal_23 <- 0 # 2 \indep 3
  n_acausal_31 <- 3 # 3 \indep 1
  
  n_no_indep <- 7 # all dependent
  
  structures <- c(
    n_all_indep,
    n_indep_1, n_indep_2, n_indep_3,
    n_causal_123, n_causal_231, n_causal_312,
    n_acausal_12, n_acausal_23, n_acausal_31,
    n_no_indep
  )
  
  structures / sum(structures)
}


CPDAG_adjacency_matrices_3var <- function() {
  all_indep <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0), 3, 3)
  
  indep_1 <- matrix(c(0, 0, 0, 0, 0, 1, 0, 1, 0), 3, 3)
  indep_2 <- matrix(c(0, 0, 1, 0, 0, 0, 1, 0, 0), 3, 3)
  indep_3 <- matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0), 3, 3)
  
  causal_123 <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), 3, 3)
  causal_231 <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
  causal_312 <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
  
  acausal_12 <- matrix(c(0, 0, 1, 0, 0, 1, 0, 0, 0), 3, 3)
  acausal_23 <- matrix(c(0, 0, 0, 1, 0, 0, 1, 0, 0), 3, 3)
  acausal_31 <- matrix(c(0, 1, 0, 0, 0, 0, 0, 1, 0), 3, 3)
  
  no_indep <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
  
  structures <- list(
    all_indep = all_indep,
    indep_1 = indep_1, indep_2 = indep_2, indep_3 = indep_3,
    causal_123 = causal_123, causal_231 = causal_231, causal_312 = causal_312, 
    acausal_12 = acausal_12, acausal_23 = acausal_23, acausal_31 = acausal_31,
    no_indep = no_indep
  )
  
  structures
}
