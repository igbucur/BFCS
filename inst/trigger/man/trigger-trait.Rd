\name{trigger.trait-methods}
\alias{trigger.trait}
\alias{trigger.trait,trigger-method}
\title{Trait-trigger analysis}
\description{
Trait-Trigger identifies, for a given trait of interest, causal gene regulator(s) that makes the trait conditionally independent of the QTL and their estimated p-values of causal regulation. These probabilities can further be used to construct a gene regulatory network.
}
\usage{
  \S4method{trigger.trait}{trigger}(triggerobj, trait, cross, thr, n.sv = NULL, addplot = TRUE)
}
\arguments{
  \item{triggerobj}{An object of class \code{\linkS4class{trigger}}. See \code{\linkS4class{trigger}} for details.}
  \item{trait}{Trait for which causal regulator is to be found. It can either be a gene-name for a gene expression trait present in \code{triggerobj} or a vector of values for the individuals present in \code{triggerobj}.}
  \item{cross}{An object of class \code{cross} obtained from \code{\link{trigger.export2cross}}. See \code{R/qtl} for more details.}
  \item{thr}{LOD threshold to search for locally linked putative causal genes (default 3).}
  \item{n.sv}{Number of surrogate variables used to model the local heterogeneity. If not set, it is computed from the expression data.}
  \item{addplot}{If \code{TRUE}, a plot of the LOD scores from a genome-scan for a single-QTL model from package \code{R/qtl}.}
}

\value{A vector of p-values associated with each tested causal regulator. 
}
\references{
Chen L.S., Emmert-Streib F., and Storey J.D. (2007) Harnessing naturally randomized transcription to infer regulatory relationships among genes. \emph{Genome Biology}, \bold{8}: R219.

Broman KW, Wu H, Sen S, Churchill GA (2003) R/qtl: QTL mapping in experimental crosses. \emph{Bioinformatics} \bold{19:} 889--890.
}
\author{
Lin S. Chen \email{lschen.stat@gmail.com}, Dipen P. Sangurdekar \email{dps@genomics.princeton.edu} and John D. Storey \email{jstorey@princeton.edu}
}
\seealso{ \code{\link{trigger.loclink}} and \code{\link{trigger.export2cross}}}
\examples{
 \dontrun{
  data(yeast)
  attach(yeast)
  triggerobj <- trigger.build(marker = marker, exp = exp, 
				   marker.pos = marker.pos, exp.pos = exp.pos)
  crossfile <- trigger.export2cross(triggerobj)
  tt.pval <- trigger.trait(triggerobj, trait = "DSE1", cross = crossfile)
  causal.reg <- names(which(p.adjust(tt.pval, method = "fdr")<.05))
  detach(yeast)
}
}
\keyword{Methods}
