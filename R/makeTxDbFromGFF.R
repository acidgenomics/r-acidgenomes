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
#' ## Ensembl
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

#' @describeIn makeTxDbFromGFF Primary function.
#' @export
makeTxDbFromGFF <- function(file, seqinfo = NULL) {
    requireNamespaces("GenomicFeatures")
    alert(sprintf(
        "Making {.var %s} from {.file %s} with {.pkg %s}::{.fun %s}.",
        "TxDb", file,
        "GenomicFeatures", "makeTxDbFromGFF"
    ))
    dataSource <- file
    if (!isAURL(dataSource)) {
        dataSource <- realpath(dataSource)
    }
    file <- .cacheIt(file)
    ## Attempt to get seqinfo automatically if not manually defined (e.g.
    ## for RefSeq). This is supported for most NCBI (e.g. GRCh38 or GRCh38.p13)
    ## and UCSC (e.g. hg38) genome builds.
    if (is.null(seqinfo)) {
        genomeBuild <- .detectGenomeBuildFromGFF(file)
        organism <- tryCatch(
            expr = detectOrganism(genomeBuild),
            error = function(e) NA_character_
        )
        seqinfo <- tryCatch(
            expr = Seqinfo(genome = genomeBuild),
            error = function(e) {
                alertWarning(paste(
                    "Automatic seqinfo detection failed.",
                    "Manual input of {.var seqinfo} is recommended."
                ))
                NULL
            }
        )
    }
    ## Additional arguments of potential future interest:
    ## - dbxrefTag: This can help override primary identifier to use.
    ## - miRBaseBuild: miRBase annotations could be useful for future genome
    ##   builds. Note that it's currently out of date with GRCh38.
    ##   https://github.com/Bioconductor/GenomicFeatures/issues/27
    args <- list(
        "file" = file,
        "dataSource" = dataSource,
        "organism" = organism
    )
    ## This step will call `.tidy_seqinfo()` to check that all chromosomes
    ## map in the seqinfo, which can be overly strict. Using a manual assignment
    ## approach below instead, inspired by similar approach in tximeta.
    ## > if (!is.null(seqinfo)) {
    ## >     args <- append(x = args, values = list("chrominfo" = seqinfo))
    ## > }
    what <- GenomicFeatures::makeTxDbFromGFF
    suppressWarnings({
        txdb <- do.call(what = what, args = args)
    })
    if (!is.null(seqinfo)) {
        assert(areIntersectingSets(names(seqinfo), names(seqinfo(txdb))))
        seqinfo(txdb) <- seqinfo[names(seqinfo(txdb))]
    }
    assert(is(txdb, "TxDb"))
    validObject(txdb)
    txdb
}
