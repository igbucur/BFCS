\name{trigger-class}
\docType{class}
\alias{trigger-class}
\title{A class to store and analyze data for Transcriptional Regulation Inference from Genetics of Gene ExpRession}
\description{
\code{trigger} is a class of objects to store and analyze data for Integrative Genomic Analysis. Use \code{\link{trigger.build}} to generate new objects of the class from input data.
}

\section{Slots}{
  \describe{
    \item{\code{exp}:}{A numeric matrix with \code{m} rows and  \code{n} columns, containing the gene expression (or intermediate trait) data.}
    \item{\code{exp.pos}:}{A matrix with \code{m} rows and 3 columns containing the chromosome number, gene start and gene end for all the genes in the gene expression matrix. The rows of \code{exp.pos} should match those of \code{exp}.}
    \item{\code{marker}:}{A matrix with \code{p} rows and  \code{n} columns, containing genotyping information.}
    \item{\code{marker.pos}:}{A matrix with \code{p} rows and 2 columns containing the chromosome number and SNP position for all the genes in the gene expression matrix. The rows of \code{exp.pos} should match those of \code{exp}.}
    \item{\code{stat}:}{A matrix of pair-wise likelihood ratio statistics for linkage analysis, with genes in rows  and markers in columns.}
    \item{\code{pvalue}:}{A matrix of parametric pvalues corresponding to statistics in the \code{stat} matrix.}
    \item{\code{mlink}:}{A list containing the results of Multi-locus linkage analysis. See \code{\link{trigger.mlink}} for details.}
    \item{\code{eqtl.R2}:}{A vector containing the proportion of genome-wide variation explained by each observed locus (eQTL).See \code{\link{trigger.eigenR2}} for details.}
    \item{\code{loc.obj}:}{A list containing the results of local-linkage probability estimation. See \code{\link{trigger.loclink}} for details.}
    }
}

\value{
An object of S4 class \code{\linkS4class{trigger}} containing the marker genotype matrix (a matrix of 1,2 for haploid genotypes, or 1,2,3 for diploid genotypes), expression matrix, marker position matrix and gene/trait position matrix with ordered coordinates in respective slots.
Use \code{slot(objectname, varname)} to retrieve individual variables from the object . Use \code{print} to see the first 10 rows and columns of the expression and marker matrix.
}


\details{
The positions in \code{marker.pos} and \code{exp.pos} matrix should be in the same units (e.g., base pair, kb, or cM). 
}
\author{
Lin S. Chen \email{lschen.stat@gmail.com}, Dipen P. Sangurdekar \email{dps@genomics.princeton.edu} and John D. Storey \email{jstorey@princeton.edu}
}
\seealso{
\code{\link{trigger.build}}, \code{\link{trigger.link}}, \code{\link{trigger.mlink}}, \code{\link{trigger.eigenR2}}, \code{\link{trigger.net}} and  \code{\link{trigger.trait}} 
}
\keyword{Classes}

