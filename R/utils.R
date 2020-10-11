#' Small internal utility functions
#'
#'
#' @author Ioan Gabriel Bucur
#' 

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

