% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R,
%   R/stripTranscriptVersions-methods.R
\name{stripTranscriptVersions}
\alias{stripTranscriptVersions}
\alias{stripTranscriptVersions,Matrix-method}
\alias{stripTranscriptVersions,character-method}
\alias{stripTranscriptVersions,matrix-method}
\title{Strip transcript identifier versions}
\usage{
stripTranscriptVersions(object, ...)

\S4method{stripTranscriptVersions}{Matrix}(object)

\S4method{stripTranscriptVersions}{character}(object)

\S4method{stripTranscriptVersions}{matrix}(object)
}
\arguments{
\item{object}{Object.}

\item{...}{Additional arguments.}
}
\value{
Modified object.
Contains transcript identifiers without version numbers.
}
\description{
Strip transcript identifier versions
}
\details{
Intentionally does not strip gene versions, if defined in the object.
To perform that action in parallel, use \code{stripVersions} instead.
}
\note{
Updated 2021-08-03.

This method is strict, and will only strip Ensembl transcript
identifiers beginning with \code{"ENS.*T"}.
}
\examples{
## Ensembl.
stripTranscriptVersions("ENST00000431238.7")

## GENCODE.
stripTranscriptVersions("ENST00000431238.7_PAR_Y")

## WormBase (no modification).
stripTranscriptVersions("cTel79B.1")
}
\seealso{
\itemize{
\item http://www.ensembl.org/info/genome/stable_ids/index.html
}
}