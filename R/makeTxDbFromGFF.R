## nolint start

#' Make TxDb from a GFF/GTF file
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
#' @param source `character(1)`.
#'   Description of data source.
#'   Defaults to file name.
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
    genomeBuild <- genome(seqinfo)[[1L]]
    if (is.null(organism)) {
        organism <- detectOrganism(genomeBuild)
    }
    assert(
        isString(organism),
        isString(genomeBuild)
    )
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
