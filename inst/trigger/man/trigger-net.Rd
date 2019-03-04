\name{trigger.net-methods}
\alias{trigger.net}
\alias{trigger.net,trigger-method}
\title{Network-Trigger analysis}
\description{
Network-Trigger analysis estimates the joint posterior probability of causal regulation for each pair of genes in the genome. These probabilities can further be used to construct a gene regulatory network.
}
\usage{
  	\S4method{trigger.net}{trigger}(triggerobj, gender = NULL, idx = NULL, 
	Bsec = 100, prob.cut = 0.7, include.loc = TRUE, seed = 123, inputfile = NULL)
}
\arguments{
  \item{triggerobj}{An object of class \code{\linkS4class{trigger}} containing slot \code{loc.obj} with local-linkage probabilities and marker indices of the best local-linkage markers for genes in the genome. See \code{\linkS4class{trigger}} and \code{\link{trigger.loclink}} for details.}
  \item{gender}{Optional. When computing statistics involving markers on sex chromosome, \code{gender} of each sample should be specified.}
  \item{idx}{Optional. One can specify the indices of selected genes as putative regulators. By default, all the genes will be selected as putative regulators.}
  \item{Bsec}{Number of iterations to perform when estimating null statistics for secondary-linkage and conditional independence.}
  \item{prob.cut}{Probability threshold. The joint regulatory probabilities  of a regulator to all the other genes will be set to zero if the local-linkage probability of the regulator is below the threshold; default \code{prob.cut = 0.7}.}
  \item{include.loc}{Logical. If \code{TRUE}, the estimated posterior probability of regulation is more conservative.}
  \item{seed}{Optional. A numeric seed for reproducible results. }
   \item{inputfile}{Optional. If provided, reads in the probability matrix from working directory.}
}
\details{
The option \code{idx} contains the indices of putative regulator genes. When the data set is large, one can use this option by selecting a subset of genes as putative regulators in one computation and parallel-computes the genome-wide regulatory probability. If \code{idx=NULL}, all the genes will be computed for probability of regulation to other genes in the data. 

If \code{include.loc = TRUE}, the joint posterior probability of regulation is the product of local-linkage, secondary-linkage and conditional independence. Otherwise, it is the product of secondary-linkage and conditional independence. The local-linkage is not a necessary condition for calculating regulation probability. If the probability of local-linkage is considered, the joint probability of regulation is more conservative. See references for details.
}
\value{A matrix of genome-wide regulatory probabilities with putative regulators in rows and regulated genes in columns. Note that the matrix is not symmetric. If gene i is estimated to be causal for gene j with high probability, the reverse is not true.
}
\references{
Chen L.S., Emmert-Streib F., and Storey J.D. (2007) Harnessing naturally randomized transcription to infer regulatory relationships among genes. \emph{Genome Biology}, \bold{8:} R219.
}
\author{
Lin S. Chen \email{lschen.stat@gmail.com}, Dipen P. Sangurdekar \email{dps@genomics.princeton.edu} and John D. Storey \email{jstorey@princeton.edu}
}
\seealso{\code{\link{trigger.loclink}}, \code{\link{trigger.netPlot2ps}} and \code{\link{trigger.trait}}}
\examples{
 \dontrun{
  data(yeast)
  attach(yeast)
  triggerobj <- trigger.build(marker = marker, exp = exp, 
            marker.pos = marker.pos, exp.pos = exp.pos)
  triggerobj <- nettrig.loc(triggerobj, window.size = 30000)
  trig.prob <- trigger.net(triggerobj, Bsec = 100)
  netPlot2ps(trig.prob)
  detach(yeast)
}
}
\keyword{Methods}
