#' Small internal utility functions
#'
#'
#' @author Ioan Gabriel Bucur
#' 

#' Derive parents list from bnlearn graph for use in compute_BGe_score.
#'  
#' @param bn_graph 
#' 
#' @return List of parents for each node, where the nodes are indexed from zero.
get_parents_list <- function(bn_graph) {
  lapply(bn_graph$nodes, function(node) which(names(bn_graph$nodes) %in% node$parents) - 1)
}
