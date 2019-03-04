\name{trigger.build}
\alias{trigger.build}
\title{Format the input data and create an Trigger object}
\description{
This function takes high-dimensional expression data and genotype data with each of their position data in the genome and creates a \code{\linkS4class{trigger}} object for subsequent analysis.
}
\usage{
  trigger.build(exp = exp,exp.pos = exp.pos, marker = marker, marker.pos = marker.pos)
}
\arguments{
  \item{exp}{A gene (or intermediate trait) by individual matrix of expression data.}
  \item{exp.pos}{A matrix containing the position information for genes (intermediate traits). The first column is the chromosome name of the gene. The second column is the starting coordinate of the gene, and the third column is the ending coordinate. Each row corresponds to one gene/trait in the \code{exp} matrix.}
  \item{marker}{A marker genotype by individual matrix.}
  \item{marker.pos}{A matrix containing the position information for markers. The first column is the chromosome name of the marker. We recommend to use integers for autosomal chromosomes and "X" for sex chromosome. The second column is the position of the marker on the chromosome. Each row corresponds to one marker in the \code{marker} matrix.}
}
\details{
The positions in \code{marker.pos} and \code{exp.pos} matrix should be in the same units (e.g., base pair, kb, or cM). 
}
\value{
An object of S4 class \code{\linkS4class{trigger}} containing the marker genotype matrix (a matrix of 1,2 for haploid genotypes, or 1,2,3 for diploid genotypes), expression matrix, marker position matrix and gene/trait position matrix with ordered coordinates in respective slots. Use \code{slot(objectname, varname)} to retrieve individual variables from the object. Use \code{\link{print}} to see the first 10 rows and columns of the expression and marker matrix.
}
\author{
Lin S. Chen \email{lschen.stat@gmail.com}, Dipen P. Sangurdekar \email{dps@genomics.princeton.edu} and John D. Storey \email{jstorey@princeton.edu}
}
\seealso{
\code{\link{trigger.link}}, \code{\link{trigger.mlink}}, \code{\link{trigger.eigenR2}}, \code{\link{trigger.net}} and  \code{\link{trigger.trait}} 
}
\examples{
\dontrun{
  data(yeast)
  attach(yeast)
  triggerobj <- trigger.build(marker = marker, exp = exp, 
                    marker.pos = marker.pos, exp.pos = exp.pos)
  print(triggerobj)
  }
}
\keyword{misc}

