#' Get RefSeq genome assembly seqinfo
#'
#' Parse the assembly report file to get `seqlengths` per chromosome.
#'
#' @export
#' @note Updated 2021-01-18.
#'
#' @param file `character(1)`.
#'   RefSeq assembly report file or URL.
#'
#' @return `Seqinfo`.
#'
#' @seealso
#' - `tximeta:::gtf2RefSeq`.
#'
#' @examples
#' file <- pasteURL(
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
#' seqinfo <- getRefSeqSeqinfo(file)
#' print(seqinfo)
getRefSeqSeqinfo <- function(file) {
    pattern <- "^([a-z0-9]+_)?GCF_[0-9]+\\.[0-9]+_(.+)_assembly_report\\.txt$"
    assert(
        isString(file),
        isMatchingRegex(pattern = pattern, x = basename(file))
    )
    file <- .cacheIt(file)
    ## e.g. GRCh38.p13, which is the format Seqinfo expects.
    ## Refer to GenomeInfoDb documentation for details on NCBI.
    genomeBuild <- sub(
        pattern = pattern,
        replacement = "\\2",
        x = basename(file)
    )
    df <- import(
        file = file,
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
