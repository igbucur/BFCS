\name{trigger.loclink-methods}
\alias{trigger.loclink}
\alias{trigger.loclink,trigger-method}
\title{Estimate local-linkage probability for each gene}
\description{
A method of class \code{\linkS4class{trigger}} to identify the best local-linkage marker for each gene and compute the local linkage probabilities.
}
\usage{
 	\S4method{trigger.loclink}{trigger}(triggerobj, gender = NULL, window.size = 30000)
}
\arguments{
  \item{triggerobj}{An object of class \code{\linkS4class{trigger}}.}
  \item{gender}{Optional. When computing linkage statistics involving markers on sex chromosome, \code{gender} of each sample should be specified.}
  \item{window.size}{Optional. The size of a window that places the putative regulator gene in the center. Every marker within the window is a candidate marker for local-linkage to the regulator gene.}
}

\value{
An updated object of class \code{\linkS4class{trigger}} containing a slot \code{loc.obj} with fields:
  \item{prob.loc}{The estimated local-linkage probability for each putative regulator gene.}
  \item{loc.idx}{The indices of the best local marker for each putative regulator gene.}

Use \code{slot(triggerobj, "loc.obj")} to retrieve the list.
}
\references{
Chen L.S., Emmert-Streib F., and Storey J.D. (2007) Harnessing naturally randomized transcription to infer regulatory relationships among genes. \emph{Genome Biology}, \bold{8:} R219.
}
\author{
Lin S. Chen \email{lschen.stat@gmail.com}, Dipen P. Sangurdekar \email{dps@genomics.princeton.edu} and John D. Storey \email{jstorey@princeton.edu}
}

\seealso{ \code{\link{trigger.trait}}}
\examples{
\dontrun{
  data(yeast)
  attach(yeast)
  triggerobj <- trigger.build(marker = marker, exp = exp, 
            marker.pos = marker.pos, exp.pos = exp.pos)
  triggerobj <- trigger.loclink(triggerobj, window.size = 30000)
  trigger.obj <- trigger.net(triggerobj, Bsec = 100)
  detach(yeast)
}
}
\keyword{Methods}

