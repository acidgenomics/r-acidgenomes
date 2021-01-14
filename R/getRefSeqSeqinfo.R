#' Get RefSeq genome assembly seqinfo
#'
#' Parse the assembly report file to get `seqlengths` per chromosome.
#'
#' @export
#' @note Updated 2021-01-14.
#'
#' @param reportFile `character(1)`.
#'   Assembly report file or URL.
#' @param genome `character(1)` or `NULL`.
#'   Genome build.
#'   If left `NULL`, will be detected from file name.
#'
#' @return `Seqinfo`.
#'
#' @seealso
#' - `tximeta:::gtf2RefSeq`.
#'
#' @examples
#' reportFile <- pasteURL(
#'     "ftp.ncbi.nlm.nih.gov",
#'     "genomes",
#'     "refseq",
#'     "vertebrate_mammalian",
#'     "Homo_sapiens",
#'     "all_assembly_versions",
#'     "GCF_000001405.38_GRCh38.p12",
#'     "GCF_000001405.38_GRCh38.p12_assembly_report.txt",
#'     protocol = "ftp"
#' )
#' seqinfo <- getRefSeqSeqinfo(reportFile = reportFile)
#' print(seqinfo)
getRefSeqSeqinfo <- function(
    reportFile,
    genomeBuild = NULL
) {
    assert(
        isString(reportFile),
        isString(genomeBuild, nullOK = TRUE)
    )
    reportFile <- .cacheIt(reportFile)
    pattern <- "^(.+)?GCF_[0-9]+\\.[0-9]+_(.+)_assembly_report.txt$"
    if (
        is.null(genomeBuild) &&
        grepl(pattern = pattern, x = basename(reportFile))
    ) {
        ## e.g. GRCh38.p13, which is the format Seqinfo expects.
        ## Refer to GenomeInfoDb documentation for details on NCBI.
        genomeBuild <- sub(
            pattern = pattern,
            replacement = "\\2",
            x = basename(reportFile)
        )
    }
    assert(isString(genomeBuild))
    df <- import(
        file = reportFile,
        format = "tsv",
        colnames = c(
            "sequenceName",
            "sequenceRole",
            "assignedMolecule",
            "assignedMoleculeLocation",
            "genbankAccn",
            "relationship",
            "refseqAccn",
            "assemblyUnit",
            "sequenceLength",
            "ucscStyleName"
        ),
        comment = "#"
    )
    cols <- c("refseqAccn", "sequenceLength")
    df <- df[, cols]
    df <- df[complete.cases(df), ]
    seqnames <- df[["refseqAccn"]]
    seqlengths <- df[["sequenceLength"]]
    assert(
        !any(is.na(seqnames)),
        !any(is.na(seqlengths)),
        hasNoDuplicates(seqnames)
    )
    seq <- Seqinfo(
        seqnames = seqnames,
        seqlengths = seqlengths,
        isCircular = NA,
        genome = genomeBuild
    )
    assert(is(seq, "Seqinfo"))
    validObject(seq)
    seq
}
