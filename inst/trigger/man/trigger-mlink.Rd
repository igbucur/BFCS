\name{trigger.mlink-methods}
\alias{trigger.mlink}
\alias{trigger.mlink,trigger-method}
\alias{Multiple-locus linkage}
\title{Multi-Locus Linkage (Epistasis) Analysis}
\description{
Multi-locus linkage (epistasis) analysis.
}
\usage{
  \S4method{trigger.mlink}{trigger}(triggerobj, prob.cut = 0.9, 
	gender = NULL, idx = NULL, B = 5, seed = 123)
}
\arguments{
  \item{triggerobj}{An object of class \code{\linkS4class{trigger}}.}
  \item{prob.cut}{Probability threshold for primary linkage.}
  \item{gender}{Optional. When computing linkage statistics involving markers on sex chromosome, \code{gender} of each sample should be specified.}
  \item{idx}{The indices for genes to be computed for multi-locus linkage.}
  \item{B}{The number of null iterations to perform.}
  \item{seed}{Optional. A numeric seed for reproducible results.}
}
\details{
When data set is large, one can the option \code{idx} to select a subset of genes in each computation and parallel-computes the genome-wide multi-locus linkage. Since the function computes the linkage probability by borrowing information across genes, at least more than 100 genes should be selected in applying this function. If \code{idx=NULL}, all the genes in the input data will be computed for multi-locus linkage.  

The current version of the function could only compute two-locus joint linkage (epistasis).  
}
\value{
An updated object of class \code{\linkS4class{trigger}} containing a slot \code{trigger.mlink} with fields:
  \item{qtl}{The major and secondary QTLs for each selected gene.}
  \item{prob}{The posterior probability of linkage for major QTL, secondary QTL, and the joint posterior probability of multi-locus linkage.} 
  \item{qvalue}{Q-value estimates for joint multi-locus linkage probabilities.}

Use \code{slot(triggerobj, "mlink")} to retrieve the list. 
}
\author{
Lin S. Chen \email{lschen.stat@gmail.com}, Dipen P. Sangurdekar \email{dps@genomics.princeton.edu} and John D. Storey \email{jstorey@princeton.edu}
}
\references{
Brem R.B., Storey J.D., Whittle J., and Kruglyak L. (2005) Genetic interactions between polymorphisms that affect gene expression in yeast. \emph{Nature}, \bold{436(7051):} 701--703.

Storey J.D., Akey J.M., and Kruglyak L. (2005) Multiple locus linkage analysis of genomewide expression in yeast. \emph{PLoS Biology}, \bold{3(8)}: 1380--1390.
}
\seealso{
 \code{\link{trigger.link} and \link{plot}} 
}
\examples{
\dontrun{
	data(yeast)
    attach(yeast)
	triggerobj <- trigger.build(marker = marker,exp = exp, 
            marker.pos = marker.pos, exp.pos = exp.pos)
	## Genome-wide multiple locus linkage analysis
	triggerobj <- trigger.mlink(triggerobj, B = 10, idx = NULL, seed = 123)
  
	plot(triggerobj, type = "trigger.mlink", qcut=0.1, bin.size=NULL)
	mlink = slot(triggerobj, "trigger.mlink")
	detach(yeast)
}
}
\keyword{Methods}

