% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapGeneNamesToNcbi.R
\name{mapGeneNamesToNcbi}
\alias{mapGeneNamesToNcbi}
\title{Map gene names to NCBI}
\usage{
mapGeneNamesToNcbi(genes, organism, taxonomicGroup = NULL, ncbi = NULL)
}
\arguments{
\item{genes}{\code{character}.
Gene names (e.g. \code{"TUT4"}).}

\item{organism}{\code{character(1)}.
Full Latin organism name (e.g. \code{"Homo sapiens"}).}

\item{taxonomicGroup}{\code{character(1)}.
NCBI FTP server taxonomic group subdirectory path (e.g. "Mammalia").
Defining this manually avoids having to query the FTP server.}

\item{ncbi}{\code{NcbiGeneIfo} or \code{NULL}.
If \code{NULL}, NCBI annotations will be downloaded automatically.}
}
\description{
Map gene names to NCBI
}
\note{
Updated 2025-04-07.
}
\examples{
## Homo sapiens.
x <- mapGeneNamesToNcbi(
    genes = c("TUT4", "ZCCHC11", "TENT3A"),
    organism = "Homo sapiens"
)
print(x)

## Mus musculus.
x <- mapGeneNamesToNcbi(
    genes = c("Nfe2l2", "Nrf2"),
    organism = "Mus musculus"
)
print(x)
}
