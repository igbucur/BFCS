% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{time_versus_sample_size}
\alias{time_versus_sample_size}
\title{A dataset containing results for reproducing Figure 9 from Section 5 in the 
IJAR article.}
\format{
A tibble data frame containing time measurements of deriving posterior
probabilities of causal regulatory relationships given sample size.
\describe{
  \item{BFCS.elapsed}{Time elapsed running BFCS}
  \item{BGe.elapsed}{Time elapsed running BGe equivalent method}
  \item{trigger.elapsed}{Time elapsed running Trigger}
  \item{size}{Network size (number of variables)}
  \item{num_edges}{Number of edges in the network}
  \item{nobs}{Number of observations generated from the network}
}
}
\usage{
time_versus_sample_size
}
\description{
A dataset containing results for reproducing Figure 9 from Section 5 in the 
IJAR article.
}
\keyword{datasets}
