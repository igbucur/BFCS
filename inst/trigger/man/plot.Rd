\name{plot}
\alias{plot}
\alias{plot,trigger,missing-method}
\title{Graphical Display of Trigger Analysis}
\description{
Graphical display of genomewide linkage map, multi-locus linkage or eQTL variation
}
\usage{
	\S4method{plot}{trigger,missing}(x,y,type = c("link", "mlink", "eqtl"),
				cutoff  = 3.3e-4, qcut = 0.1, bin.size = NULL)
}
\arguments{
	\item{x}{An object of class \code{\linkS4class{trigger}}.}
  \item{y}{Ignore option, not used.}
 	\item{type}{An argument describing the type of plot. Select from \code{link} (default) for genome-wide linkage map, \code{eqtl.R2} for graphical display of eQTL-R^2 contribution or \code{mlink} for display of genome-wide epistasis effect.}
  	\item{cutoff}{Threshold value for \code{link}. The measures below the threshold are called significant and are plotted. }
  	\item{qcut}{Q-value threshold for \code{mlink}. The joint multi-locus linkage probabilities with q-values below the threshold are called significant and are plotted. }
  	\item{bin.size}{Optional for \code{mlink}. If not \code{NULL}, each chromosome will be divided into several bins, each with size \code{bin.size}. Markers within a bin will be considered as at a same position. }
}

\author{
Lin S. Chen \email{lschen.stat@gmail.com}, Dipen P. Sangurdekar \email{dps@genomics.princeton.edu} and John D. Storey \email{jstorey@princeton.edu}
}
\seealso{
 	\code{\link{trigger.link}}, \code{\link{trigger.mlink}} and \code{\link{trigger.eigenR2}} 
}
\examples{
 \dontrun{
  
  data(yeast)
  attach(yeast)
  triggerobj <- trigger.build(marker=marker, exp=exp, 
				marker.pos=marker.pos, exp.pos=exp.pos)
  triggerobj <- trigger.link(triggerobj, gender=NULL, norm=TRUE)
  plot(triggerobj,type = "link", cutoff=1e-5)
  triggerobj <- trigger.eigenR2(triggerobj, adjust=FALSE)
  plot(triggerobj, type = "eigenR2")
  triggerobj<- trigger.mlink(triggerobj, B=5, seed=123)
  plot(triggerobj, qcut=0.1, bin.size=NULL)
  
  detach(yeast)

}
}
\keyword{Methods}

