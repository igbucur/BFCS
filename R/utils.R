#' Small internal utility functions and global variables
#' @author Ioan Gabriel Bucur
IJAR_theme <- theme_classic() + theme(
  text = element_text(size = 15), 
  legend.text = element_text(size = 15), 
  legend.position = "bottom",
  legend.key.size = unit(1.5, "lines"), 
  legend.key.width = unit(1.5, "lines"),
  aspect.ratio = 1, plot.margin = unit(c(0, 0, 0, 0), "lines"))

#' Derive parents list from bnlearn graph for use in compute_BGe_score.
#'  
#' @param bn_graph object of class "bn" representing a bnlearn graph
#' 
#' @return List of parents for each node, where the nodes are indexed from zero.
get_parents_list <- function(bn_graph) {
  lapply(bn_graph$nodes, function(node) which(names(bn_graph$nodes) %in% node$parents) - 1)
}


#' Function for extracting off-diagonal elements of a matrix.
#'
#' @param mat Matrix of any type.
#'
#' @return Vector containing off-diagonal elements of matrix.
#' 
#' @keywords internal
off_diagonal <- function(mat) {
  mat[lower.tri(mat) | upper.tri(mat)]
}


#' This function returns a list of all 11 different Markov equivalence classes
#' for causal graphs (DAGs) on three variables.
#'
#' @param format The format in which to return the list. Currently, only "amat"
#' (adjacency matrix) and "bnlearn" (class "bn") formats are supported.
#'
#' @return A list containing the 11 different Markov equivalence classes.
#' @export
#'
#' @examples 
#' get_Markov_equivalence_class_list(format = 'bnlearn') # default
#' get_Markov_equivalence_class_list(format = 'amat')
get_Markov_equivalence_class_list <- function(format = 'bnlearn') {
  
  # DAGs on triples (bnlearn format)
  
  if (format == 'bnlearn') {
    
    # 1. Empty (1 DAG w/ BK)
    empty <- bnlearn::empty.graph(c("L_k", "T_i", "T_j"))
    
    # 2. Independent (4 DAG w/ BK)
    indep_1_ij <- bnlearn::set.arc(empty, "T_i", "T_j")
    indep_1_ji <- bnlearn::set.arc(empty, "T_j", "T_i")
    indep_2 <- bnlearn::set.arc(empty, "L_k", "T_j")
    indep_3 <- bnlearn::set.arc(empty, "L_k", "T_i")
    
    # 3. Causal (3 DAG w/ BK)
    causal_123 <- bnlearn::set.arc(indep_2, "L_k", "T_i")
    causal_231 <- bnlearn::set.arc(indep_3, "T_i", "T_j")
    causal_312 <- bnlearn::set.arc(indep_2, "T_j", "T_i")
    
    # 4. Acausal (2 DAG w / BK)
    acausal_12 <- bnlearn::set.arc(indep_2, "T_i", "T_j")
    acausal_23 <- bnlearn::set.arc(bnlearn::set.arc(empty, "T_i", "L_k"), "T_j", "L_k")
    acausal_31 <- bnlearn::set.arc(indep_3, "T_j", "T_i")
    
    # 5. Full (2 DAG w / BK)
    full_ij <- bnlearn::set.arc(acausal_12, "L_k", "T_i")
    full_ji <- bnlearn::set.arc(acausal_31, "L_k", "T_j")
    
    mec_list <- list(
      empty,
      indep_1_ij,
      indep_2,
      indep_3,
      causal_123,
      causal_231,
      causal_312,
      acausal_12,
      acausal_23,
      acausal_31,
      full_ij
    )
  } else if(format == 'amat') {
    
    # DAGs on triples (amat format, upper diagonal)
    
    # 1. Empty (1 DAG w/ BK)
    empty <- matrix(0, 3, 3)
    
    # 2. Independent (4 DAG w/ BK)
    indep_1 <- empty; indep_1[2, 3] <- 1 # alternatively indep_1[3, 2] <- 1
    indep_2 <- empty; indep_2[1, 3] <- 1
    indep_3 <- empty; indep_3[1, 2] <- 1
    
    # 3. Causal (3 DAG w/ BK)
    causal_123 <- indep_2; causal_123[1, 2] <- 1
    causal_231 <- indep_3; causal_231[2, 3] <- 1
    causal_312 <- indep_2; causal_312[3, 2] <- 1
    
    # 4. Acausal (2 DAG w / BK)
    acausal_12 <- indep_2; acausal_12[2, 3] <- 1
    acausal_23 <- empty; acausal_23[2, 1] <- acausal_23[3, 1] <- 1
    acausal_31 <- indep_3; acausal_31[3, 2] <- 1
    
    # 5. Full (2 DAG w / BK)
    full <- acausal_12; full[1, 2] <- 1
    # alternatively full <- acausal_31; full[1, 3] <- 1 
    
    # Markov equivalence classes with BK
    mec_list <- list(
      empty, 
      indep_1,
      indep_2,
      indep_3,
      causal_123,
      causal_231,
      causal_312,
      acausal_12,
      acausal_23,
      acausal_31,
      full
    )
  } else {
    stop('Unknown format. Currently, only "bnlearn" and "amat" formats are
         supported.')
  }
  
  names(mec_list) <- c(
    'empty',
    'indep_1', 'indep_2', 'indep_3',
    'causal_123', 'causal_231', 'causal_312',
    'acausal_12', 'acausal_23', 'acausal_31',
    'full'
  )
  
  mec_list
}

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
