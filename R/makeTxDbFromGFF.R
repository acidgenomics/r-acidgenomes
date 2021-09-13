## NOTE This currently isn't working as desired for UCSC GTF files.
##      TxDb doesn't keep track of transcript identifiers correctly.
##      Should contain: ENST00000450305



## nolint start

#' Make TxDb from a GFF/GTF file
#'
#' Wrapper for GenomicFeatures `makeTxDbFromGFF` importer.
#'
#' @name makeTxDbFromGFF
#' @note Updated 2021-08-09.
#'
#' @inheritParams AcidRoxygen::params
#' @inheritParams params
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
#' ## GENCODE ====
#' ## > gtfFile <- pasteURL(
#' ## >     "ftp.ebi.ac.uk",
#' ## >     "pub",
#' ## >     "databases",
#' ## >     "gencode",
#' ## >     "Gencode_human",
#' ## >     "release_36",
#' ## >     "gencode.v36.annotation.gtf.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > txdb <- AcidGenomes::makeTxDbFromGFF(file = gtfFile)
#' ## > print(txdb)
#'
#' ## RefSeq ====
#' ## > gffFile <- pasteURL(
#' ## >     "ftp.ncbi.nlm.nih.gov",
#' ## >     "genomes",
#' ## >     "refseq",
#' ## >     "vertebrate_mammalian",
#' ## >     "Homo_sapiens",
#' ## >     "all_assembly_versions",
#' ## >     "GCF_000001405.38_GRCh38.p12",
#' ## >     "GCF_000001405.38_GRCh38.p12_genomic.gff.gz",
#' ## >     protocol = "ftp"
#' ## > )
#' ## > txdb <- AcidGenomes::makeTxDbFromGFF(file = gffFile)
#' ## > print(txdb)
NULL

## nolint end

#' @describeIn makeTxDbFromGFF Primary function.
#' @export
makeTxDbFromGFF <- function(file) {
    pkgs <- .packages()
    assert(.isSupportedGFF(file))
    ## Check for input of unsupported files.
    ## See `.gffPatterns` for details.
    denylist <- c(
        "refseq_gtf" = paste0(
            "^([0-9a-z]_)?",
            "(GC[AF]_[0-9]+\\.[0-9]+)",
            "_([^_]+)",
            "_(.+)",
            "\\.gtf",
            "(\\.gz)?$"
        )
    )
    if (isMatchingRegex(
        pattern = denylist[["refseq_gtf"]],
        x = basename(file)
    )) {
        abort(sprintf(
            paste(
                "Unsupported file: {.file %s}.",
                "Use RefSeq GFF instead of GTF.",
                sep = "\n"
            ),
            basename(file)
        ))
    }
    alert(sprintf(
        "Making {.cls %s} from {.file %s} with {.pkg %s}::{.fun %s}.",
        "TxDb", file,
        "GenomicFeatures", "makeTxDbFromGFF"
    ))
    requireNamespaces("GenomicFeatures")
    if (isAFile(file)) {
        file <- realpath(file)
    }
    meta <- getGFFMetadata(file)
    organism <- meta[["organism"]]
    seqinfo <- .getSeqinfo(meta)
    assert(
        isString(organism),
        isAny(seqinfo, c("Seqinfo", "NULL"))
    )
    ## Additional arguments of potential future interest:
    ## - dbxrefTag: This can help override primary identifier to use.
    ## - miRBaseBuild: miRBase annotations could be useful for future genome
    ##   builds. Note that it's currently out of date with GRCh38.
    ##   https://github.com/Bioconductor/GenomicFeatures/issues/27
    args <- list(
        "file" = .cacheIt(file),
        "dataSource" = file,
        "organism" = organism
    )
    if (!is.null(seqinfo)) {
        args <- append(x = args, values = list("chrominfo" = seqinfo))
    }
    what <- GenomicFeatures::makeTxDbFromGFF
    suppressWarnings({
        txdb <- do.call(what = what, args = args)
    })
    assert(is(txdb, "TxDb"))
    ## Stash the GFF metadata, so we can access in `makeGRangesFromGFF()`.
    attr(txdb, which = "gffMetadata") <- meta
    validObject(txdb)
    forceDetach(keep = pkgs)
    txdb
}
