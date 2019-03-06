#' Vectorized compute Bayes factors of covariance structures using the Bayesian
#' Gaussian equivalent (BGe) score.
#'
#' @param suffstat_list List of sufficient statistics for a data set, containing
#' the means of the variables and their covariance matrix.
#' @param no_samples number of samples
#' @param graph_list List containing objects of class "bn", each of which
#' represents a graph structure for which we compute the BGe score.
#'
#' @return Matrix of size (number datasets) x (number graph structures) that
#' contains the BGe score for each combination of (dataset, graph).
#' @export
#'
#' @examples 
#' # 1. Some idealized examples
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
#' compute_BGe_score_vectorized(suffstat_list, no_samples = 10000, mec_list)
#' 
#' # 2. Example on yeast data set
#' library(trigger)
#' data(yeast)
#' 
#' L <- yeast$marker
#' Trait <- yeast$exp
#' yeast_data <- t(rbind(yeast$marker, yeast$exp))
#' means.LT <- colMeans(yeast_data)
#' cov.LT <- cov(yeast_data)
#' 
#' # select triples 
#' my_grid <- subset(expand.grid(L_k = 1:10, T_i = 11:20, T_j = 11:20), T_i != T_j)
#' 
#' suffstat_list <- apply(my_grid, 1, function(row) {
#'   list(means = means.LT[row], covmat = cov.LT[row, row])
#' })
#'
# compute_BGe_score_vectorized(suffstat_list, 112, get_Markov_equivalence_class_list())
compute_BGe_score_vectorized <- function(
  suffstat_list,
  no_samples,
  graph_list
) {
  
  # Rcpp:::sourceCpp('src/compute_BGe_score.cpp')
  
  lev <- matrix(0, length(suffstat_list), length(graph_list))

  for (j in 1:length(graph_list)) {
    
    graph <- graph_list[[j]]
    
    # TODO: find a more elegant solution
    # provide BGe_score with a list of parents for all nodes
    parents_list <- lapply(graph$nodes, function(node) which(names(graph$nodes) %in% node$parents) - 1)
    
    for (i in 1:length(suffstat_list)) {
      
      ss <- suffstat_list[[i]]
      
      lev[i, j] <- compute_BGe_score(no_samples, ss$means, ss$covmat, 
                        parents_list, alpha_w = 5, nu = rep(0, 3))
    }
  }
  
  lev
}


# simpler local score function for root nodes.
bge_root_node_vectorized <- function(
  grid, means, covmat, no_samples, var, alpha_mu = 1,
  nu = structure(rep(0, ncol(grid)), names = ncol(grid)), alpha_w = ncol(grid) + 2) {
  
  N <- no_samples
  n <- ncol(grid)
  
  stopifnot(var %in% 1:n)
  
  # first term.
  res <- rep(0.5 * (log(alpha_mu) - log(N + alpha_mu)), nrow(grid))
  
  # Gamma_l ratio in the second term.
  res <- res + lgamma((N + alpha_w - n + 1)/2) - lgamma((alpha_w - n + 1)/2)
  
  # leftover from the second term.
  res = res - N/2 * log(pi)
  
  # third term, numerator.
  t = (alpha_mu + alpha_w - n - 1)/(alpha_mu + 1)
  res = res + (alpha_w - n + 1)/2 * log(t)
  
  # third term, denominator.
  nu <- as.numeric(nu[var])
  
  logR <- apply(grid, 1, function(row) {
    #select <- c(node, parents)
    #local = data[, c(node, parents), drop = FALSE]
    R <- t + covmat[row, row][var, var] * (N - 1) +
      (N * alpha_w) / (N + alpha_w) * (means[row][var] - nu)^2
    
    log(R)
  })
  
  res = res - (N + alpha_w - n + 1)/2 * logR
  
  return(res)
  
}#BGE_ROOT_NODE

# NOTE: compute_BGe_score_vectorized is still faster
# my_grid <- subset(expand.grid(L_k = 1:10, T_i = 11:20, T_j = 11:20), T_i != T_j)
# 
# microbenchmark(
# res1 <- bge_root_node_vectorized(my_grid, means.LT, cov.LT, 112, 1) +
#   bge_root_node_vectorized(my_grid, means.LT, cov.LT, 112, 2) +
#   bge_root_node_vectorized(my_grid, means.LT, cov.LT, 112, 3),
#   compute_BGe_score_vectorized(suffstat_list, 112, mec_list[1]))


bge_conditional_node_vectorized <- function(
  grid, means, covmat, no_samples, node, parents, alpha_mu = 1,
  nu = structure(rep(0, ncol(grid)), names = ncol(grid)), alpha_w = ncol(grid) + 2) {
  
  N <- no_samples
  n <- ncol(grid)
  p <- length(parents)
  
  stopifnot(node %in% 1:n)
  stopifnot(all(parents %in% 1:n))
  
  # first term.
  res <- rep(0.5 * (log(alpha_mu) - log(N + alpha_mu)), nrow(grid))
  
  # Gamma ratio in the second term.
  res <- res + lgamma((N + alpha_w - n + p + 1)/2) - lgamma((alpha_w - n + p + 1)/2)
  
  # leftover from the second term.
  res <- res - N/2 * log(pi)
  
  # third term, ratio of the determinants of the prior T matrices.
  t <- (alpha_mu + alpha_w - n - 1)/(alpha_mu + 1)
  
  res = res + (alpha_w - n + p + 1)/2 * (p + 1) * log(t) -
    (alpha_w - n + p)/2 * (p) * log(t)
  
  select <- c(node, parents)
  
  # third term, ratio of the determinants of the posterior R matrices.
  nu <- nu[select]
  
  logR <- apply(grid, 1, function(row) {
    #select <- c(node, parents)
    #local = data[, c(node, parents), drop = FALSE]
    R <- diag(t, p + 1) + covmat[row, row][select, select] * (N - 1) +
      (N * alpha_w) / (N + alpha_w) *
      outer(means[row][select] - nu, means[row][select] - nu)
    
    c(log(det(R)), log(det(R[-1, -1])))
  })
  
  
  
  res = res + (N + alpha_w - n + p)/2 * logR[2,] -
    (N + alpha_w - n + p + 1)/2 * logR[1, ]
  
  return(res)
  
}#BGE_CONDITIONAL_NODE



