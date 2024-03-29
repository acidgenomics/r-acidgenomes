% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downloadUcscGenome.R
\name{downloadUcscGenome}
\alias{downloadUcscGenome}
\title{Download UCSC reference genome}
\usage{
downloadUcscGenome(
  organism,
  genomeBuild = NULL,
  outputDir = getwd(),
  cache = FALSE
)
}
\arguments{
\item{organism}{\code{character(1)}.
Full Latin organism name (e.g. \code{"Homo sapiens"}).}

\item{genomeBuild}{\code{character(1)}.
UCSC genome build assembly name (e.g. \code{"hg38"}).
If set \code{NULL}, defauls to the most recent build available.}

\item{outputDir}{\code{character(1)}.
Output directory path.}

\item{cache}{\code{logical(1)}.
Cache URLs locally, using BiocFileCache internally.}
}
\value{
Invisible \code{list}.
}
\description{
Download UCSC reference genome
}
\note{
Updated 2023-11-22.
}
\section{Genome}{

\itemize{
\item \verb{<GENOME_BUILD>.chrom.sizes}: Two-column tab-separated text file containing
assembly sequence names and sizes.
\item \verb{<GENOME_BUILD>.chromAlias.txt}: Sequence name alias file, one line for
each sequence name.  First column is sequence name followed by tab
separated alias names.
}
}

\section{Transcriptome}{

\itemize{
\item \code{mrna.fa.gz}: Human mRNA from GenBank. This sequence data is updated
regularly via automatic GenBank updates.
\item \code{refMrna.fa.gz}: RefSeq mRNA from the same species as the genome.
This sequence data is updated regularly via automatic GenBank updates.
}
}

\section{Gene annotations}{


This directory contains GTF files for the main gene transcript sets where
available. They are sourced from the following gene model tables:
\code{knownGene} (GENCODE) and \code{ncbiRefSeq} (NCBI RefSeq).
}

\examples{
## This example is bandwidth intensive.
## > downloadUcscGenome(organism = "Homo sapiens")
}
\seealso{
\itemize{
\item \href{http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/}{hg38 bigZips}.
\item \href{https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/}{hg38 pinned analysis set (for NGS pipelines)}.
\item \href{https://genome.ucsc.edu/FAQ/FAQformat.html}{UCSC genome FAQ}.
\item \href{https://genome.ucsc.edu/FAQ/FAQdownloads.html}{UCSC download FAQ}.
\item \href{http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format}{UCSC genes in GTF or GFF format wiki}.
}
}
