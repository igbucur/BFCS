#' Compute the posterior probabilities of L_k -> T_i -> T_j for a simulated GRN 
#' using bfcs::compute_BGe_score_vectorized.
#'
#' @param GRN List produced by simulate_GRN containing GRN structure and statistics.
#' @param prior Prior on causal graphs (default uniform on DMAGs with bkg. info.)
#'
#' @return A numeric matrix of posterior probabilities for each possible causal
#' relationship between the expression traits. 
#' @export
#'
#' @examples
#' GRN <- simulate_GRN(10, 10, 100, 0.1)
#' get_BGe_probabilities(GRN)
get_BGe_probabilities <- function(GRN, prior = uniform_prior_GRN_DMAG()) {
  
  ngen <- nrow(GRN$L)
  nexp <- nrow(GRN$Trait)
  
  # Markov equivalence classes for three nodes
  mec_list <- get_Markov_equivalence_class_list()
  
  bnlearn_sim <- lapply(1:nexp, function(i) {
    my_grid <- expand.grid(
      L_k = 1:ngen,
      T_i = ngen + i,
      T_j = (ngen + 1):(ngen + nexp)
    )
    
    suffstat_list <- apply(my_grid, 1, function(row) {
      list(means = GRN$means[row], covmat = GRN$cov.matrix[row, row])
    }) 
    
    compute_BGe_score_vectorized(suffstat_list, ncol(GRN$L), mec_list)
  })
  
  sapply(1:nexp, function(i) {
    Bf <- apply(bnlearn_sim[[i]], 1, function(lev) exp(lev - max(lev)) * prior)
    prob <- matrix(Bf[6, ] / colSums(Bf), ngen, nexp)
    prob <- apply(prob, 2, max)
    prob[i] <- 0
    prob
  })
}

#' Compute the posterior probabilities of L_k -> T_i -> T_j for a simulated GRN
#' using the BiDAG package.
#'
#' @param GRN List produced by simulate_GRN containing GRN structure and statistics.
#' @param prior Prior on causal graphs (default uniform on DMAGs with bkg. info.)
#'
#' @return A numeric matrix of posterior probabilities for each possible causal
#' relationship between the expression traits. 
#' @export
#'
#' @examples
#' GRN <- simulate_GRN(5, 5, 100, 0.1)
#' get_BiDAG_probabilities(GRN)
get_BiDAG_probabilities <- function(GRN, prior = uniform_prior_GRN_DMAG()) {
  
  ngen <- nrow(GRN$L)
  nexp <- nrow(GRN$Trait)
  
  prob_matrix <- matrix(0, nexp, nexp)
  
  mec_list <- get_Markov_equivalence_class_list(format = "amat")
  
  for (i in 1:nexp) {
    T_i <- GRN$Trait[i, ]
    
    for (j in 1:nexp) {
      if (i == j) next
      T_j <- GRN$Trait[j, ]
      
      probs <- sapply(1:ngen, function(k) {
        
        L_k <- GRN$L[k, ]
        
        triple <- data.frame(L_k, T_i, T_j)
        # triple <- apply(triple, 2, scale) # normalize
        
        # Compute posterior probabilities
        lev <- sapply(mec_list, function(mec) {
          BiDAG::DAGscore(3, BiDAG::scoreparameters(3, "bge", triple), mec)
        })
        
        ev <- exp(lev - max(lev)) * prior
        ev <- ev / sum(ev) # normalize probabilities
        ev[6] # probability of L_k -> T_i -> T_j
      })
      prob_matrix[j, i] <- max(probs)
    }
  }
  
  prob_matrix
}