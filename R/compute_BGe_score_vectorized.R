#' Vectorized compute Bayes factors of covariance structures using the Bayesian
#' Gaussian equivalent (BGe) score.
#'
#' @param suffstat_list List of sufficient statistics for a data set, containing
#' the means of the variables and their covariance matrix.
#' @param num_samples number of samples
#' @param graph_list List containing objects of class "bn", each of which
#' represents a graph structure for which we compute the BGe score.
#'
#' @return Matrix of size (number datasets) x (number graph structures) that
#' contains the BGe score for each combination of (dataset, graph).
#' @export
#'
#' @examples 
#' mec_list <- get_Markov_equivalence_class_list(format = "bnlearn")
#' suffstat_list = list(
#'   # three independent variables, zero mean
#'   list( 
#'     means = rep(0, 3),
#'     covmat = diag(3)
#'   ),
#'   # three independent variables, mean one
#'   list(
#'     means = rep(1, 3),
#'     covmat = diag(3)
#'   ),
#'   # a single variable repeated three times (correlation one)
#'   list(
#'     means = rep(0, 3),
#'     covmat = matrix(1, 3, 3)
#'   )
#' )
#' compute_BGe_score_vectorized(suffstat_list, num_samples = 10000, mec_list)
compute_BGe_score_vectorized <- function(
  suffstat_list,
  num_samples,
  graph_list
) {
  
  # Rcpp:::sourceCpp('src/compute_BGe_score.cpp')
  
  log_evidence <- matrix(0, length(suffstat_list), length(graph_list))

  for (j in 1:length(graph_list)) {
    
    graph <- graph_list[[j]]
    
    # TODO: find a more elegant solution
    # provide BGe_score with a list of parents for all nodes
    parents_list <- lapply(graph$nodes, function(node) which(names(graph$nodes) %in% node$parents) - 1)
    
    for (i in 1:length(suffstat_list)) {
      
      ss <- suffstat_list[[i]]
      
      log_evidence[i, j] <- compute_BGe_score(num_samples, ss$means, ss$covmat, 
                        parents_list, alpha_w = 5, nu_vec = rep(0, 3))
    }
  }
  
  log_evidence
}


