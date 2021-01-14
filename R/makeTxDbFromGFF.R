## nolint start

#' Make TxDb from GRanges
#'
#' @name makeTxDbFromGFF
#' @note Updated 2021-01-14.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @param seqinfo `Seqinfo`.
#'   Information about the chromosomes.
#' @param source `character(1)`.
#'   Description of data source.
#'   Defaults to file name.
#'
#' @details
#' This step can be noisy and generate expected warnings, which are suppressed:
#'
#' ```
#' some exons are linked to transcripts not found in the file
#' The following orphan exon were dropped
#' The following orphan CDS were dropped
#' ```
#'
#' @return `TxDb`.
#'
#' @seealso
#' - `GenomicFeatures::makeTxDbFromGFF()`.
#' - `GenomicFeatures::supportedMiRBaseBuildValues()`.
#'   Note that *Homo sapiens* GRCh38 isn't currently supported in mirbase.db.
#' - [TxDb.Hsapiens.UCSC.hg38.knownGene](https://bioconductor.org/packages/TxDb.Hsapiens.UCSC.hg38.knownGene/).
#' - https://stackoverflow.com/questions/38603668
#' - https://stackoverflow.com/questions/16517795
#' - https://github.com/Bioconductor/GenomicFeatures/blob/master/R/makeTxDbFromGRanges.R
#' - https://github.com/Bioconductor/GenomicFeatures/issues/26
#'
#' @examples
#' file <- pasteURL(AcidGenomesTestsURL, "ensembl.gtf")
#' txdb <- makeTxDbFromGFF(file)
#' print(txdb)
NULL

## nolint end



#' @describeIn makeTxDbFromGFF Primary function.
#' @export
makeTxDbFromGFF <- function(
    file,
    seqinfo,  ## Pass this in for RefSeq.
    organism = NULL,
    genomeBuild = NULL,
    source = file
) {
    requireNamespaces("GenomicFeatures")
    assert(
        isString(file),
        is(seqinfo, "Seqinfo"),
        isString(organism, nullOK = TRUE),
        isString(genomeBuild, nullOK = TRUE),
        isString(source)
    )
    alert(sprintf(
        "Making {.var %s} from {.file %s} with {.pkg %s}::{.fun %s}.",
        "TxDb", file,
        "GenomicFeatures", "makeTxDbFromGFF"
    ))

    ## FIXME ATTEMPT TO DETECT ORGANISM AND GENOME BUILD AUTOMATICALLY FROM
    ## FILE NAME.

    assert(
        isString(organism),
        isString(genomeBuild)
    )
    file <- .cacheIt(file)
    args <- list(
        "file" = file,
        ## Information about the chromosomes.
        "chrominfo" = seqinfo,
        ## Description string about the origin of the data file.
        ## Can use URL input here.
        "dataSource" = source,
        ## Handle either GFF or GTF input automatically.
        "format" = "auto"

        ## > dataSource = NA,
        ## > organism = "FIXME",
        ## > circ_seqs = NULL,
        ## > chrominfo = "FIXME"
        ## This could be useful for future genome builds with miRs.
        ## > miRBaseBuild = NA,
        ## > metadata = NULL,
        ## This can help override feature name ("e.g. GeneID") used.
        ## FIXME DO WE NEED TO RETHINK THIS FOR TRANSCRIPT IDs?
        ## > dbxrefTag = "GeneID"
    )
    what <- GenomicFeatures::makeTxDbFromGFF
    suppressWarnings({
        txdb <- do.call(what = what, args = args)
    })
    assert(is(txdb, "TxDb"))
    validObject(txdb)
    txdb
}



#' @describeIn makeTxDbFromGFF Alias for GTF files.
#' @export
makeTxDbFromGTF <- makeTxDbFromGFF
