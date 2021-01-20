## FIXME ATTEMPT TO DETECT THE GENOME BUILD FROM THE FILE HEADER...



## nolint start

#' Make TxDb from a GFF/GTF file
#'
#' Wrapper for GenomicFeatures `makeTxDbFromGFF` importer.
#'
#' @name makeTxDbFromGFF
#' @note Updated 2021-01-20.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @param seqinfo `Seqinfo` or `NULL`.
#'   **Recommended.** Information about the chromosomes.
#'
#' @return `TxDb`.
#'
#' @seealso
#' - `GenomicFeatures::makeTxDbFromGFF()`.
#' - `GenomicFeatures::supportedMiRBaseBuildValues()`.
#'   Note that *Homo sapiens* GRCh38 isn't currently supported in mirbase.db.
#' - [TxDb.Hsapiens.UCSC.hg38.knownGene](https://bioconductor.org/packages/TxDb.Hsapiens.UCSC.hg38.knownGene/).
#'
#' @examples
#' ## UCSC.
#'
#'
#' ## RefSeq.
#' gffFile <- pasteURL(
#'     "ftp.ncbi.nlm.nih.gov",
#'     "genomes",
#'     "refseq",
#'     "vertebrate_mammalian",
#'     "Homo_sapiens",
#'     "all_assembly_versions",
#'     "GCF_000001405.38_GRCh38.p12",
#'     "GCF_000001405.38_GRCh38.p12_genomic.gff.gz",
#'     protocol = "ftp"
#' )
#' reportFile <- sub(
#'     pattern = "_genomic\\.gff\\.gz$",
#'     replacement = "_assembly_report.txt",
#'     x = gffFile
#' )
#' seqinfo <- getRefSeqSeqinfo(reportFile)
#' txdb <- makeTxDbFromGFF(file = gffFile, seqinfo = seqinfo)
#' print(txdb)
NULL

## nolint end



## FIXME FOR SEQINFO INPUT, DO WE NEED TO MATCH SEQLEVELS?

#' @describeIn makeTxDbFromGFF Primary function.
#' @export
makeTxDbFromGFF <- function(file, seqinfo = NULL) {
    requireNamespaces("GenomicFeatures")
    assert(
        isString(file),
        is(seqinfo, "Seqinfo") || is.null(seqinfo)
    )
    alert(sprintf(
        "Making {.var %s} from {.file %s} with {.pkg %s}::{.fun %s}.",
        "TxDb", file,
        "GenomicFeatures", "makeTxDbFromGFF"
    ))
    dataSource <- file
    if (!isAURL(dataSource)) {
        dataSource <- realpath(dataSource)
    }
    if (is.null(seqinfo)) {
        alertWarning(paste(
            "Input of {.var seqinfo} is recommended.",
            "This helps define chromosome seqlengths and genome metadata."
        ))
        genomeBuild <- NA_character_
        organism <- NA_character_
    } else {
        ## e.g. "GRCh38.p12".
        genomeBuild <- genome(seqinfo)[[1L]]
        if (is.null(organism)) {
            ## e.g. "Homo sapiens".
            organism <- tryCatch(
                expr = detectOrganism(genomeBuild),
                error = function(e) NA_character_
            )
        }
    }
    file <- .cacheIt(file)
    ## Additional arguments of potential future interest:
    ## - dbxrefTag: This can help override primary identifier to use.
    ## - miRBaseBuild: miRBase annotations could be useful for future genome
    ##   builds. Note that it's currently out of date with GRCh38.
    ##   https://github.com/Bioconductor/GenomicFeatures/issues/27
    args <- list(
        "file" = file,
        "dataSource" = dataSource,
        "format" = "auto",
        "organism" = organism
    )
    if (!is.null(seqinfo)) {
        args <- append(
            x = args,
            values = list(
                "chrominfo" = seqinfo,
                "circ_seqs" = isCircular(seqinfo),
            )
        )
    }
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
