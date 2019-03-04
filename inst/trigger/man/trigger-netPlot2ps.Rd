\name{trigger.netPlot2ps-methods}
\alias{trigger.netPlot2ps}
\alias{trigger.netPlot2ps,trigger-method}
\title{Write the network from a trigger probability matrix to a postscript file}
\description{
Write the network from a trigger probability matrix to a postscript file.
}
\usage{
  \S4method{trigger.netPlot2ps}{trigger}(triggerobj, trig.prob, filenam = NULL, pcut = 0.95, 
	layout = c("radial", "energy-minimized", "circular","hierarchical"),
	 node.color = NULL, edge.color = NULL, node.shape = NULL, nreg = 20)
}
\arguments{
  \item{triggerobj}{An object of class \code{\linkS4class{trigger}}.}
  \item{trig.prob}{A network-Trigger regulatory probability matrix with putative regulator genes in rows and putative regulated genes in columns. See \code{\link{trigger.net}} for details.}
  \item{filenam}{The output file name, without extension. If the name is not specified, the network will be write to the files \code{temp.ps} and \code{temp.dot} at the current directory.}
  \item{pcut}{Threshold value for regulatory probabilities. The probabilities above the threshold are called significant and the corresponding regulatory relationships are plotted.}
  \item{layout}{The layout of the output network. One can choose from "\code{radial}" (default), "\code{energy-minimized}", "\code{circular}" or "\code{hierarchical}" layouts. You can specify just the initial letter.}
  \item{node.color}{The color of the nodes (genes). The default color is green.}
  \item{edge.color}{The color of the edges. The default color is blue.}
  \item{node.shape}{The shape of nodes (genes) if the number of regulatory relationships is below 1000. If that number is above 1000, the shape of nodes will be \code{dot}.}
  \item{nreg}{The number of top regulators to be selected. These selected top regulators will be plotted in red ellipses with their gene names labeled inside.}
}
\details{
To use this function, please install the software \code{Graphviz}, which is available at http://www.graphviz.org/. For large networks, layout "\code{radial}" or "\code{energy-minimized}" is recommended. If the total number of significant regulatory relationships (directed edges) of the network is below 1000, we plot each node (gene) as a "\code{box}" with its name labeled inside. Otherwise, we plot each gene as a "\code{dot}" without name labeled to facilitate visualization. The top \code{nreg} (by default \code{nreg = 20}) regulators will be plotted in red ellipses labeled with their names. 

See manual of \code{Graphviz} for other available colors and shapes of nodes.
}
\author{
Lin S. Chen \email{lschen.stat@gmail.com}, Dipen P. Sangurdekar \email{dps@genomics.princeton.edu} and John D. Storey \email{jstorey@princeton.edu}
}
\seealso{
 \code{\link{trigger.link} and \link{trigger.mlink}} 
}
\examples{
\dontrun{
 data(yeast)
  attach(yeast)
  triggerobj <- trigger.build(marker = marker, exp = exp, 
				marker.pos = marker.pos, exp.pos = exp.pos)
  triggerobj <- trigger.loclink(triggerobj, window.size = 30000)
  trig.prob <- trigger.net(triggerobj, Bsec = 100)

  trigger.netPlot2ps(trig.prob, pcut = 0.95, layout = "e", filenam = "net95", nreg = 20)
  
  detach(yeast)
}
}
\keyword{misc}

