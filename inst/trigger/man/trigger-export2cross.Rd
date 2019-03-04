\name{trigger.export2cross-methods}
\alias{trigger.export2cross}
\alias{trigger.export2cross,trigger-method}
\title{Export Trigger data to R/qtl's cross class object}
\description{
\code{trigger.export2cross} exports \code{\linkS4class{trigger}} data from \code{triggerobj} to a \code{cross} format for Trait-Trigger analysis. See \code{\link{trigger.trait}} for details.
}
\usage{
  \S4method{trigger.export2cross}{trigger}(triggerobj, plotarg = TRUE, verbose = TRUE, warning = FALSE)
}
\arguments{
  \item{triggerobj}{An object of class \code{\linkS4class{trigger}}.}
  \item{plotarg}{Logical. If \code{TRUE}, the function plots the default plot from the \code{R/qtl} package while reading in the genotype data.}
  \item{verbose}{Logical. If \code{TRUE}, the function lists the default output from the \code{R/qtl} package while reading in the genotype data.}
  \item{warning}{Logical. If \code{FALSE}, the function suppresses warnings output from the \code{R/qtl} package while reading in the genotype data.}
}
\details{
The \code{trigger.export2cross} command writes a \code{csv} format file \dQuote{geno_trait_data.csv} to the working directory and reads it using the \code{read.cross} command.
}
\value{An object of class cross from the \code{R/qtl} package.
}
\references{
Broman KW, Wu H, Sen S, Churchill GA (2003) R/qtl: QTL mapping in experimental crosses. \emph{Bioinformatics} \bold{19:} 889--890. 
}
\author{
Lin S. Chen \email{lschen.stat@gmail.com}, Dipen P. Sangurdekar \email{dps@genomics.princeton.edu} and John D. Storey \email{jstorey@princeton.edu}
}
\seealso{ \code{\link{trigger.trait}}}
\examples{
  data(yeast)
  attach(yeast)
  triggerobj <- trigger.build(marker = marker, exp = exp, 
					marker.pos = marker.pos, exp.pos = exp.pos)
  crossfile <- trigger.export2cross(triggerobj, plotarg = TRUE, verbose = TRUE, warning = FALSE)
  tt.pval <- trigger.trait(triggerobj, trait  =  "DSE1", cross  =  crossfile)
  causal.reg <- names(which(p.adjust(tt.pval, method = "fdr")<.05))
  detach(yeast)

}
\keyword{Methods}
