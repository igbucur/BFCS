\name{trigger.eigenR2-methods}
\alias{trigger.eigenR2}
\alias{trigger.eigenR2,trigger-method}
\title{Estimate the proportion of genome-wide variation explained by each eQTL}
\description{
Estimate eqtl-R2, the proportion of genome-wide variation explained by each eQTL and identify linkage hotspots.
}
\usage{
	\S4method{trigger.eigenR2}{trigger}(triggerobj, adjust = FALSE, meanR2 = FALSE)
}
\arguments{
  \item{triggerobj}{An object of class \code{\linkS4class{trigger}}.}
  \item{adjust}{Logical. If \code{TRUE}, the estimated R-square for each locus will be adjusted for small sample size effect. Recommend to use when sample size is less than 100.}
  \item{meanR2}{Logical. If \code{TRUE}, the function computes the mean of R-squares of genome-wide gene expression for each locus.}
}

\value{
An updated object of class \code{\linkS4class{trigger}} with a slot \code{loc.obj} containing the proportion of genome-wide variation explained by each observed locus (eQTL). Use \code{slot(triggerobj, "eigenR2")} to retrieve the eqtl-R2 values as a vector.
}

\author{
Lin S. Chen \email{lschen.stat@gmail.com}, Dipen P. Sangurdekar \email{dps@genomics.princeton.edu} and John D. Storey \email{jstorey@princeton.edu}
}
\references{
Chen L.S. and Storey J.D. (2008) Eigen-R2 for dissecting variation in high-dimensional studies. \emph{Bioinformatics} \bold{24(19):} 2260--2262.
}
\seealso{
 \code{\link{plot}} 
}
\examples{
\dontrun{
  data(yeast)
  attach(yeast)
  triggerobj <- trigger.build(marker = marker, exp = exp, 
			marker.pos = marker.pos, exp.pos = exp.pos)
  triggerobj <- trigger.eigenR2(triggerobj, adjust = FALSE)
  plot(triggerobj, type = "eigenR2")
  eqtlR2 <- slot(triggerobj, "eigenR2")
  detach(yeast)

}
}
\keyword{Methods}

