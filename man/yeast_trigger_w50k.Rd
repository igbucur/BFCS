% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{yeast_trigger_w50k}
\alias{yeast_trigger_w50k}
\title{A dataset containing trigger posterior probabilities on the yeast dataset,
which we ran using a windows size of 50000 for trigger.loclink,
for reproducing Table 2 in Section 4 of the IJAR article, as well as Tables 1 
and 2 from Section 4 in the PGM 2018 article.}
\format{
6216x6216 matrix, where the element [i, j] indicates the estimated
probability that T_i -> T_j
}
\usage{
yeast_trigger_w50k
}
\description{
A dataset containing trigger posterior probabilities on the yeast dataset,
which we ran using a windows size of 50000 for trigger.loclink,
for reproducing Table 2 in Section 4 of the IJAR article, as well as Tables 1 
and 2 from Section 4 in the PGM 2018 article.
}
\keyword{datasets}
