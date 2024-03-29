% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/export-methods.R
\name{export}
\alias{export}
\alias{export,TxToGene,character-method}
\title{Export data from R}
\usage{
export(object, con, ...)

\S4method{export}{TxToGene,character}(object, con, ...)
}
\arguments{
\item{object}{Object.}

\item{con}{\code{character(1)} or \code{connection}.
Data connection. Use \code{character(1)} to represent a file path.}

\item{...}{Additional arguments.}
}
\value{
Invisible \code{character} or \code{list} of file paths.
}
\description{
Export data from R
}
\details{
The \code{TxToGene} method automatically disables writing of column names, which
is the intended input format for tximport.
}
\note{
Updated 2023-09-20.
}
\examples{
object <- S4Vectors::DataFrame(
    "txId" = c(
        "tx0001",
        "tx0002",
        "tx0003",
        "tx0004"
    ),
    "geneId" = c(
        "gene0001",
        "gene0001",
        "gene0002",
        "gene0002"
    )
)
object <- TxToGene(object)
con <- file.path(AcidBase::tempdir2(), "tx2gene.csv")
export(object = object, con = con)
x <- readLines(con, n = 4L)
print(x)
AcidBase::unlink2(con)
}
