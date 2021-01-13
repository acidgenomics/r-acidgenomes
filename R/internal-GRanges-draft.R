#' Get RefSeq seqinfo
#'
#' Parse the assembly report file to get seqlengths per chromosome.
#'
#' @note Updated 2021-01-13.
#' @noRd
#'
#' @rdname reportFile `character(1)`.
#'   Assembly report file or URL.
#' @rdname genome `character(1)` or `NULL`.
#'   Genome build.
#'   If left `NULL`, will be detected from file name.
#'
#' @seealso
#' - `tximeta:::gtf2RefSeq`.
#'
#' @examples
#' reportFile <- pasteURL(
#'     "ftp.ncbi.nlm.nih.gov",
#'     "genomes",
#'     "all",
#'     "GCF",
#'     "000",
#'     "001",
#'     "405",
#'     "GCF_000001405.39_GRCh38.p13",
#'     "GCF_000001405.39_GRCh38.p13_assembly_report.txt",
#'      protocol = "ftp"
#' )
#' object <- .getRefSeqSeqinfo(reportFile = reportFile)
.getRefSeqSeqinfo <- function(
    reportFile,
    genomeBuild = NULL
) {
    assert(isString(reportFile))
    pattern <- "^GCF_[0-9]+\\.[0-9]+_(.+)_assembly_report.txt$"
    if (
        is.null(genomeBuild) &&
        grepl(pattern = pattern, x = basename(reportFile))
    ) {
        ## e.g. GRCh38.p13, which is the format Seqinfo expects.
        ## Refer to GenomeInfoDb documentation for details on NCBI.
        genomeBuild <- sub(
            pattern = pattern,
            replacement = "\\1",
            x = basename(reportFile)
        )
    }
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
