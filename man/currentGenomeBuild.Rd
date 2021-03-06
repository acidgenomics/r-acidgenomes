% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/currentGenomeBuild.R
\name{currentGenomeBuild}
\alias{currentGenomeBuild}
\alias{currentEnsemblGenomeBuild}
\alias{currentGencodeGenomeBuild}
\alias{currentRefSeqGenomeBuild}
\alias{currentUCSCGenomeBuild}
\title{Current genome build}
\usage{
currentEnsemblGenomeBuild(organism)

currentGencodeGenomeBuild(organism)

currentRefSeqGenomeBuild(organism, taxonomicGroup = NULL)

currentUCSCGenomeBuild(organism)
}
\arguments{
\item{organism}{\code{character(1)}.
Full Latin organism name (e.g. \code{"Homo sapiens"}).}

\item{taxonomicGroup}{\code{character(1)}.
\emph{Only applies to RefSeq}.
FTP server taxonomic group subdirectory path (e.g. "vertebrate_mammalian").
Defining this manually avoids having to query the FTP server.}
}
\value{
\code{character(1)}.
Genome assembly build version.
}
\description{
Fetch the current genome build (assembly) version from online resources.
}
\note{
Updated 2021-04-27.
}
\examples{
currentEnsemblGenomeBuild("Homo sapiens")
currentGencodeGenomeBuild("Homo sapiens")
currentRefSeqGenomeBuild(
    organism = "Homo sapiens",
    taxonomicGroup = "vertebrate_mammalian"
)
currentUCSCGenomeBuild("Homo sapiens")
}
\seealso{
\itemize{
\item \href{https://rest.ensembl.org}{Ensembl REST API}.
\item \href{https://genome.ucsc.edu/goldenPath/help/api.html}{UCSC REST API}.
\item \href{https://ftp.ncbi.nlm.nih.gov/genomes/refseq/}{RefSeq genomes}.
}
}
