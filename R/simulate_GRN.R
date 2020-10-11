#' Generate gene regulatory network (GRN) data from structural equation model.
#'
#' @param ngen Integer number of genetic markers.
#' @param nexp Integer number of expression traits.
#' @param nobs Integer number of observations in the generated data set.
#' @param prob_edge Numeric probability of having a causal interaction (edge)
#' between two expression traits. Passed on to \link[pcalg]{randomDAG}.
#' @param seed Integer seed for reproducible random number generation.
#' @param lwr_edge_weight Lower limit of edge weight (causal strength T_i -> T_j)
#' @param upr_edge_weight Upper limit of edge weight (causal strength T_i -> T_j)
#' 
#' @return A data set consisting of \emph{nobs} samples from a random generated 
#' gene regulatory network. This includes \emph{ngen} binomially distributed 
#' variables representing the genetic markers associated with the genetic makeup 
#' of each individual observation. It also includes \emph{nexp} columns 
#' representing the gene expression levels of \emph{nexp} traits, which have an 
#' underlying causal graph describing their interactions.
#' @export
#'
#' @examples simulate_GRN(10, 10, 100, 0.1, 1634)
simulate_GRN <- function(ngen, nexp, nobs, prob_edge,  seed = NULL,
                         lwr_edge_weight = -1, upr_edge_weight = 1) {
  
  set.seed(seed)
  
  # Generate causal graph expressing causal regulatory relationships
  G <- pcalg::randomDAG(nexp, prob_edge, V = paste0("T", 1:nexp),
                        lB = lwr_edge_weight, uB = upr_edge_weight)
  
  # Generate random minor allele frequencies (MAF)
  theta <- stats::runif(ngen, 0.1, 0.5) 
  
  # Generate genetic marker data from a binomial distribution
  L <- t(sapply(theta, stats::rbinom, n = nobs, size = 1)) + 1
  rownames(L) <- paste0("L", 1:ngen)
  # we consider a single chromosome, locus identical to row number
  L.pos <- matrix(1, ngen, 2); L.pos[, 2] <- 1:ngen 
  
  B <- t(methods::as(G, "matrix")) 
  print(paste("Number of edges:", sum(B != 0)))
  C <- matrix(0, nexp, ngen); C[1:ngen, 1:ngen] <- diag(ngen)
  
  # Generate expression data on nexp traits
  Trait <- solve(diag(nexp) - B) %*% (C %*% L + matrix(stats::rnorm(nexp * nobs), nexp, nobs))
  T.pos <- matrix(1, nexp, 3); T.pos[, 2:3] <- 1:nexp
  
  # combine marker and trait data to compute summary statistics
  data <- t(rbind(L, Trait))
  
  list(L = L, Trait = Trait, L.pos = L.pos, T.pos = T.pos, 
       cov.matrix = stats::cov(data), means = colMeans(data),
       cor.matrix = stats::cor(data), B = B, num_edges = length(G@edgeData@data))
}
