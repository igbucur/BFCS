#' A dataset containing results for reproducing Figure 3 in Section 4 of the 
#' PGM 2018 and IJAR articles.
#'
#' @format A three-level list containing the Bayes factors obtained for the six 
#' scenarios described in Figures 3a-f. Each scenario has a list of nine matrices
#' with Bayes factors derived on increasing number of samples.
"Bayes_factors_consistency"

#' A dataset containing results for reproducing Figure 4, as well as Tables 2 
#' and 3 from Section 4 in the PGM 2018 article.
#'
#' @format A three-level list containing the computed posterior probabilities
#' for two different simulated GRNs (sparse, dense) from which (100, 1000) 
#' samples were generated, obtained using four different algorithms 
#' (trigger, BFCS DAG, BFCS DMAG, BGe).
"simulated_GRN_probabilities_PGM"

#' A dataset containing results for reproducing Figures 4-7, as well as Table 2
#' from Section 4 in the IJAR article.
#'
#' @format A three-level list containing the computed posterior probabilities
#' for three different simulated GRNs (sparse, less dense, more dense) from
#' which (100, 1000) samples were generated, obtained using five different 
#' algorithms (trigger, BFCS DAG, BFCS DMAG, BFCS loclink, BGe).
"simulated_GRN_probabilities_IJAR"

#' A dataset containing BFCS DMAG posterior probabilities on the yeast dataset
#' for reproducing Table 2 in Section 4 of the IJAR article, as well as Tables 1 
#' and 2 from Section 4 in the PGM 2018 article.
#'
#' @format 6216x6216 matrix, where the element [i, j] indicates the estimated
#' probability that T_i -> T_j
"yeast_BFCS_DMAG"

#' A dataset containing posterior probabilities on the yeast dataset,
#' as reported in the supplement of [Chen et al. (2008)]
#' (https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-10-r219),
#' for reproducing Table 2 in Section 4 of the IJAR article, as well as Tables 1 
#' and 2 from Section 4 in the PGM 2018 article.
#'
#' @format 6216x6216 matrix, where the element [i, j] indicates the estimated
#' probability that T_i -> T_j
"yeast_trigger_Chen"

#' A dataset containing trigger posterior probabilities on the yeast dataset,
#' which we ran using a windows size of 50000 for trigger.loclink,
#' for reproducing Table 2 in Section 4 of the IJAR article, as well as Tables 1 
#' and 2 from Section 4 in the PGM 2018 article.
#'
#' @format 6216x6216 matrix, where the element [i, j] indicates the estimated
#' probability that T_i -> T_j
"yeast_trigger_w50k"

#' A dataset containing results for reproducing Figure 8 from Section 5 in the 
#' IJAR article.
#'
#' @format A tibble data frame containing time measurements of deriving posterior
#' probabilities of causal regulatory relationships given network size.
#' \describe{
#'   \item{BFCS.elapsed}{Time elapsed running BFCS}
#'   \item{BGe.elapsed}{Time elapsed running BGe equivalent method}
#'   \item{trigger.elapsed}{Time elapsed running Trigger}
#'   \item{size}{Network size (number of variables)}
#'   \item{num_edges}{Number of edges in the network}
#' }
"time_versus_network_size"

#' A dataset containing results for reproducing Figure 9 from Section 5 in the 
#' IJAR article.
#'
#' @format A tibble data frame containing time measurements of deriving posterior
#' probabilities of causal regulatory relationships given sample size.
#' \describe{
#'   \item{BFCS.elapsed}{Time elapsed running BFCS}
#'   \item{BGe.elapsed}{Time elapsed running BGe equivalent method}
#'   \item{trigger.elapsed}{Time elapsed running Trigger}
#'   \item{size}{Network size (number of variables)}
#'   \item{num_edges}{Number of edges in the network}
#'   \item{nobs}{Number of observations generated from the network}
#' }
"time_versus_sample_size"

