% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{uniform_prior_GRN_DMAG}
\alias{uniform_prior_GRN_DMAG}
\title{Uniform prior on DMAG structures with background knowledge.}
\usage{
uniform_prior_GRN_DMAG()
}
\value{
Logical vector indicating the probability of Markov Equivalence
Classes over three variables, when X_1 precedes (X_2, X_3), assuming a
uniform prior over DMAG causal structures.
}
\description{
Uniform prior on DMAG structures with background knowledge.
}
\examples{
prior <- uniform_prior_GRN_DMAG()
}
