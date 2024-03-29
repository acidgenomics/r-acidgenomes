% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emptyRanges.R
\name{emptyRanges}
\alias{emptyRanges}
\title{Generate empty genomic ranges}
\usage{
emptyRanges(
  names,
  seqname = c("unknown", "transgene", "spike"),
  mcolnames = NULL
)
}
\arguments{
\item{names}{\code{character}.
Gene or transcript names.}

\item{seqname}{\code{character(1)}.
Name of the alternative chromosome to be defined in \code{seqnames} where these
ranges will be grouped. Defaults to \code{"unknown"} but \code{"transgene"}
(transgenes) and \code{"spike"} (spike-ins) are also supported.}

\item{mcolnames}{\code{character} or \code{NULL}.
Metadata column names to be defined in the \code{mcols} of the \code{GRanges} return.
Normally this does not need to be defined; useful when combining with another
\code{GRanges} that contains metadata.}
}
\value{
\code{GRanges}.
}
\description{
Utility function that provides support for creating internal \code{GRanges} for
transgene and FASTA spike-in sequences.
}
\note{
Updated 2019-08-21.
}
\examples{
## Unknown/dead genes.
emptyRanges("ENSG00000000000", seqname = "unknown")

## Transgenes.
emptyRanges(c("EGFP", "TDTOMATO", "GAL4"), seqname = "transgene")

## Spike-ins.
emptyRanges("ERCC", seqname = "spike")
}
\seealso{
\code{help("seqinfo", "GenomeInfoDb")}.
}
