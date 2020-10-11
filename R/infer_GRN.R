#' Compute the posterior probabilities of L_k -> T_i -> T_j for a simulated GRN 
#' using BFCS::compute_Bayes_factors_vectorized.
#'
#' @param GRN List produced by simulate_GRN containing GRN structure and statistics.
#' @param prior Prior on causal graphs (default uniform on DMAGs with bkg. info.)
#' @param loc_link Integer vector specifying genetic marker with strongest local
#' linkage for each expression trait. If not specified, we use the default BFCS 
#' maximum probability search strategy.
#'
#' @return A numeric matrix of posterior probabilities for each possible causal
#' relationship between the expression traits. 
#' @export
#'
#' @examples
#' GRN <- simulate_GRN(10, 10, 1000, 0.1)
#' compute_BFCS_probabilities(GRN)
#' compute_BFCS_probabilities(GRN, loc_link = 10:1)
compute_BFCS_probabilities <- function(GRN, prior = uniform_prior_GRN_DMAG(), loc_link = NULL) {
  
  ngen <- nrow(GRN$L)
  nexp <- nrow(GRN$Trait)
  
  # For fairness in calculating time complexity, we recompute GRN$cor.matrix
  GRN$cor.matrix <- stats::cor(t(rbind(GRN$L, GRN$Trait)))
  
  cLT <- GRN$cor.matrix[1:ngen, (ngen+1):(ngen+nexp)]
  cTT <- GRN$cor.matrix[(ngen+1):(ngen+nexp), (ngen+1):(ngen+nexp)]
  
  Bayes_factors <- lapply(1:nexp, function(i) {
    compute_Bayes_factors_vectorized(rep(cLT[, i], nexp), c(cLT), rep(cTT[i, ], each = ngen),
                                     num_samples = ncol(GRN$L))
  })
  
  posterior_probabilities <- sapply(1:nexp, function(i) {
    Bf <- prior * t(Bayes_factors[[i]])
    prob <- matrix(Bf[6, ] / colSums(Bf), ngen, nexp)
    
    if (is.null(loc_link)) {
      prob <- apply(prob, 2, max)
    } else {
      prob <- prob[loc_link[i], ]
    }

    prob[i] <- 0
    prob
  })
  
  dimnames(posterior_probabilities) <- list(
    rownames(GRN$Trait), rownames(GRN$Trait)
  )
  
  posterior_probabilities
}


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
#' compute_BGe_probabilities(GRN)
compute_BGe_probabilities <- function(GRN, prior = uniform_prior_GRN_DMAG()) {
  
  ngen <- nrow(GRN$L)
  nexp <- nrow(GRN$Trait)
  
  # Markov equivalence classes for three nodes
  mec_list <- get_Markov_equivalence_class_list()
  
  Bayes_factors <- lapply(1:nexp, function(i) {
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
  
  posterior_probabilities <- sapply(1:nexp, function(i) {
    Bf <- apply(Bayes_factors[[i]], 1, function(lev) exp(lev - max(lev)) * prior)
    prob <- matrix(Bf[6, ] / colSums(Bf), ngen, nexp)
    prob <- apply(prob, 2, max)
    prob[i] <- 0
    prob
  })
  
  dimnames(posterior_probabilities) <- list(
    rownames(GRN$Trait), rownames(GRN$Trait)
  )
  
  posterior_probabilities
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
#' compute_BiDAG_probabilities(GRN)
compute_BiDAG_probabilities <- function(GRN, prior = uniform_prior_GRN_DMAG()) {
  
  ngen <- nrow(GRN$L)
  nexp <- nrow(GRN$Trait)
  
  posterior_probabilities <- matrix(0, nexp, nexp)
  
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
          BiDAG::DAGscore(BiDAG::scoreparameters(3, "bge", triple), mec)
        })
        
        ev <- exp(lev - max(lev)) * prior
        ev <- ev / sum(ev) # normalize probabilities
        ev[6] # probability of L_k -> T_i -> T_j
      })
      posterior_probabilities[j, i] <- max(probs)
    }
  }
  
  dimnames(posterior_probabilities) <- list(
    rownames(GRN$Trait), rownames(GRN$Trait)
  )
  
  posterior_probabilities
}

#' Compute the posterior probabilities of L_k -> T_i -> T_j for a simulated GRN
#' using the trigger package.
#'
#' @param GRN List produced by simulate_GRN containing GRN structure and statistics.
#' @param trigger_obj Precomputed object of class 'trigger' (build + loclink).
#' @param seed Random seed for reproducibility.
#' @param window.size Window size when looking for genetic markers. Default
#' option is value chosen in Chen et al. (2007)
#'
#' @return A numeric matrix of posterior probabilities for each possible causal
#' relationship between the expression traits. 
#' @export
#'
#' @examples
#' GRN <- simulate_GRN(100, 100, 100, 0.01)
#' compute_trigger_probabilities(GRN)
compute_trigger_probabilities <- function(GRN, trigger_obj = NULL, seed = NULL, window.size = 50000) {
  
  if (is.null(trigger_obj)) {
    trigger_obj <- trigger::trigger.build(marker = GRN$L, exp = GRN$Trait, 
                                 marker.pos = GRN$L.pos, exp.pos = GRN$T.pos)
    
    utils::capture.output(
      trigger_obj <- trigger::trigger.loclink(trigger_obj, window.size = window.size)
    )
  }
  
  tryCatch({
    oldwd <- setwd("~/")
    t(trigger::trigger.net(trigger_obj, seed = seed))
  }, error = function(e) {
    print(e)
    NULL
  }, finally = {
    # Clean up file produced when running trigger.net
    suppressWarnings(file.remove('~/net_trigg_prob.txt'))
    setwd(oldwd)
  })
}
