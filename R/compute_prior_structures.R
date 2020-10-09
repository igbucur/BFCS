#' Uniform prior on DAG structures with background knowledge.
#'
#' @return Logical vector indicating the probability of Markov Equivalence
#' Classes over three variables, when X_1 precedes (X_2, X_3), assuming a
#' uniform prior over DAG causal structures.
#' @export
#'
#' @examples
#' prior <- uniform_prior_GRN_DAG()
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

#' Uniform prior on DMAG structures with background knowledge.
#'
#' @return Logical vector indicating the probability of Markov Equivalence
#' Classes over three variables, when X_1 precedes (X_2, X_3), assuming a
#' uniform prior over DMAG causal structures.
#' @export
#'
#' @examples
#' prior <- uniform_prior_GRN_DMAG()
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