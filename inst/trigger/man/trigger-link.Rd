\name{trigger.link-methods}
\alias{trigger.link}
\alias{trigger.link,trigger-method}
\title{Genomewide eQTL analysis}
\description{
A method of class \code{\linkS4class{trigger}} for genomewide Expression-trait QTL analysis. This function estimates the linkage statistic and parametric p-value for each gene expression to every locus in the genome. 
}
\usage{
  \S4method{trigger.link}{trigger}(triggerobj, gender = NULL, norm = TRUE)
}
\arguments{
  \item{triggerobj}{An object of class \code{\linkS4class{trigger}}.}
  \item{gender}{Optional. When computing linkage statistics involving markers on sex chromosome, \code{gender} of each sample should be specified.}
  \item{norm}{Logical. If \code{TRUE}, each row of expression matrix \code{exp} in the \code{triggerobj} will be transformed to follow a standard normal distribution, based on the rank of value.}
}
\value{
An updated object of class \code{\linkS4class{trigger}} containing slots:
  \item{stat}{A matrix of pair-wise likelihood ratio statistics for linkage analysis, with genes in rows and markers in columns.}
  \item{pvalue}{A matrix of parametric pvalues corresponding to statistics in the \code{stat} matrix.}
Use \code{slot(triggerobj, "stat")} and \code{slot(triggerobj, "pvalue")} to retrieve the values.
}
\author{
Lin S. Chen \email{lschen.stat@gmail.com}, Dipen P. Sangurdekar \email{dps@genomics.princeton.edu} and John D. Storey \email{jstorey@princeton.edu}
}
\seealso{
 \code{\link{plot} and \link{trigger.mlink}} 
}
\examples{
\dontrun{
  data(yeast)
  attach(yeast)
  triggerobj <- trigger.build(marker = marker, exp = exp, 
				marker.pos = marker.pos, exp.pos = exp.pos)
  triggerobj <- trigger.link(triggerobj, gender = NULL, norm = TRUE)
  plot(triggerobj,type = "link", cutoff = 1e-5)
  stat = slot(triggerobj, "stat"); pvalue = slot(triggerobj, "pvalue")
  detach(yeast)

}
}
\keyword{Methods}
