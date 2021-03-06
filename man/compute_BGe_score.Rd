% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{compute_BGe_score}
\alias{compute_BGe_score}
\title{Compute Bayesian Gaussian equivalent (BGe) score for a given Bayesian network
on given sufficient statistics from data.}
\usage{
compute_BGe_score(N, means, covmat, parents, alpha_w, nu_vec, alpha_mu = 1)
}
\arguments{
\item{N}{Integer number of observations.}

\item{means}{Numeric vector containing variable means.}

\item{covmat}{Numeric covariance matrix.}

\item{parents}{List of parents for each node to describe Bayesian network.}

\item{alpha_w}{Numeric parameter of the covariance's Wishart prior (scale matrix).}

\item{nu_vec}{Numeric parameter of the covariance's Wishart prior (degrees of freedom)}

\item{alpha_mu}{Numeric parameter of the mean's Gaussian prior (expected value).}
}
\value{
Bayesian Gaussian equivalent (BGe) score for given data on the
given Bayesian network.
}
\description{
Compute Bayesian Gaussian equivalent (BGe) score for a given Bayesian network
on given sufficient statistics from data.
}
\details{
For a better understanding of how the BGe score is obtained see the article 
[Geiger and Heckerman (2002)](https://projecteuclid.org/euclid.aos/1035844981)
as well as the discussion in [Kuipers et al. (2014)](https://projecteuclid.org/euclid.aos/1407420013).
}
