% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EnsemblGenes.R
\name{EnsemblGenes}
\alias{EnsemblGenes}
\title{Generate an \code{EnsemblGenes} object}
\usage{
EnsemblGenes(object)
}
\arguments{
\item{object}{\code{GRanges}.
Must contain Ensembl annotations previously generated by
\code{makeGRangesFromEnsembl} or \code{makeGRangesFromGff}.}
}
\description{
Generate an \code{EnsemblGenes} object
}
\note{
Updated 2023-04-26.
}
\examples{
## > object <- makeGRangesFromEnsembl("Homo sapiens")
## > object <- as(object, "GRanges")
## > object <- EnsemblGenes(object)
## > print(object)
}
