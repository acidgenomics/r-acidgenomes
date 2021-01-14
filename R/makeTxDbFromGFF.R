## nolint start

#' Make TxDb from a GFF/GTF file
#'
#' Wrapper for GenomicFeatures `makeTxDbFromGFF` importer.
#'
#' @name makeTxDbFromGFF
#' @note Updated 2021-01-14.
#' @note For Ensembl and GENCODE genomes, consider using
#'   `makeEnsDbFromGFF` instead.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
#'
#' @param seqinfo `Seqinfo`.
#'   Information about the chromosomes.
#' @param source `character(1)` or `NULL`.
#'   Description of data source.
#'   If left `NULL`, defaults to file name.
#'
#' @details
#' This step can be noisy and generate expected warnings, which are
#' intentionally suppressed:
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
#'
#' @examples
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



#' @describeIn makeTxDbFromGFF Primary function.
#' @export
makeTxDbFromGFF <- function(
    file,
    seqinfo,
    organism = NULL,
    source = NULL
) {
    requireNamespaces("GenomicFeatures")
    assert(
        isString(file),
        is(seqinfo, "Seqinfo"),
        isString(organism, nullOK = TRUE),
        isString(genomeBuild, nullOK = TRUE),
        isString(source, nullOK = TRUE)
    )
    alert(sprintf(
        "Making {.var %s} from {.file %s} with {.pkg %s}::{.fun %s}.",
        "TxDb", file,
        "GenomicFeatures", "makeTxDbFromGFF"
    ))
    ## e.g. "GRCh38.p12".
    genomeBuild <- genome(seqinfo)[[1L]]
    if (is.null(organism)) {
        ## e.g. "Homo sapiens".
        organism <- detectOrganism(genomeBuild)
    }
    assert(
        isString(organism),
        isString(genomeBuild)
    )
    if (is.null(source)) {
        ## Ensure we resolve the file path, if defining as data source.
        if (!isAURL(file)) file <- realpath(file)
        source <- file
    }
    file <- .cacheIt(file)
    args <- list(
        "file" = file,
        "chrominfo" = seqinfo,
        "dataSource" = source,
        "format" = "auto",
        "organism" = organism
        ## Don't think this is necessary when passing in seqinfo.
        ## > "circ_seqs" = isCircular(seqinfo),
        ## This can help override feature name ("e.g. GeneID") used.
        ## > dbxrefTag = "GeneID"
        ## miRBase annotations could be useful for future genome builds.
        ## Note that it's currently out of date with GRCh38.
        ## https://github.com/Bioconductor/GenomicFeatures/issues/27
        ## > miRBaseBuild = NA
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
