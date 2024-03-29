% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/TxToGene-methods.R
\name{TxToGene}
\alias{TxToGene}
\alias{TxToGene,DFrame-method}
\alias{TxToGene,GRanges-method}
\alias{TxToGene,GRangesList-method}
\alias{TxToGene,data.frame-method}
\alias{TxToGene,matrix-method}
\title{Transcript-to-gene mappings}
\usage{
TxToGene(object, ...)

\S4method{TxToGene}{DFrame}(object, quiet = FALSE)

\S4method{TxToGene}{GRanges}(object, ...)

\S4method{TxToGene}{GRangesList}(object, ...)

\S4method{TxToGene}{data.frame}(object, ...)

\S4method{TxToGene}{matrix}(object, ...)
}
\arguments{
\item{object}{Object.}

\item{quiet}{\code{logical(1)}.
Perform command quietly, suppressing messages.}

\item{...}{Arguments pass through to \code{DFrame} method.}
}
\value{
\code{TxToGene}.
}
\description{
Transcript-to-gene mappings
}
\note{
Updated 2023-11-29.
}
\examples{
## DFrame ====
df <- S4Vectors::DataFrame(
    "txId" = c(
        "ENST00000635602.1",
        "ENST00000635506.1"
    ),
    "geneId" = c(
        "ENSG00000283061.1",
        "ENSG00000283061.1"
    )
)
object <- TxToGene(df)
print(object)
}
\seealso{
\code{\link[=makeTxToGene]{makeTxToGene()}}.
}
