#' Generate gene regulatory network (GRN) data from structural equation model.
#'
#' @param ngen Integer number of genetic markers.
#' @param nexp Integer number of expression traits.
#' @param nobs Integer number of observations in the generated data set.
#' @param prob_edge Numeric probability of having a causal interaction (edge)
#' between two expression traits. Passed on to \link[pcalg]{randomDAG}.
#' @param seed Integer seed for random number generation.
#'
#' @return A data set consisting of nobs samples from a random generated gene
#' regulatory network. This includes \emph{ngen} binomially distributed variables
#' representing the genetic markers associated with the genetic makeup of each
#' individual observation. It also includes \emph{nexp} columns representing 
#' the gene expression levels of \emph{nexp} traits, which have an underlying 
#' causal graph describing their interactions.
#' @export
#'
#' @examples simulate_GRN(10, 10, 100, 0.1, 1634)
simulate_GRN <- function(ngen, nexp, nobs, prob_edge, seed = NULL) {
  
  set.seed(seed)
  
  # dense graph
  G <- pcalg:::randomDAG(nexp, prob_edge, lB = 1, uB = 1, V = paste0("T", 1:nexp))
  
  stopifnot(nexp >= ngen)
  
  ## Generate data
  
  theta <- runif(ngen, 0.1, 0.5) 
  
  L <- t(sapply(theta, rbinom, n = nobs, size = 1)) + 1
  rownames(L) <- paste0("L", 1:ngen)
  L.pos <- matrix(1, ngen, 2); L.pos[, 2] <- 1:ngen # one chromosome, each in the same position
  
  
  B <- t(as(G, "matrix")); print(paste("Number of edges:", sum(B)))
  C <- matrix(0, nexp, ngen); C[1:ngen, 1:ngen] <- diag(ngen)
  
  Tr <- solve(diag(nexp) - B) %*% (C %*% L + matrix(rnorm(nexp * nobs), nexp, nobs))
  # Tr <- t(apply(Tr, 1, function(row) row / sd(row))) # normalize to Var[X] = 1
  
  Trait <- Tr
  Trait <- t(apply(Tr, 1, function(x) qnorm( rank(x) / (length(x) + 1) ) ) ) # re-normalize data
  T.pos <- matrix(1, nexp, 3); T.pos[, 2:3] <- 1:nexp
  
  data <- t(rbind(L, Tr))
  cor.matrix <- cor(data) # unnormalized
  
  list(L = L, Trait = Trait, L.pos = L.pos, T.pos = T.pos, 
       cor.matrix = cor.matrix, B = B)
}