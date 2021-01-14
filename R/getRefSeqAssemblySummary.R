## nolint start

#' Get the RefSeq assembly metadata
#'
#' @section Stable reference assembly:
#' The latest reference assembly, linked under the `reference/`
#' subdirectory changes over time and is not considered "stable". For improved
#' reproduciblity, track a reference one version behind
#' (e.g. use "GCF_000001405.39_GRCh38.p12" instead of current
#' "GCF_000001405.39_GRCh38.p13" build).
#'
#' This approach to maintaining a current reference build differs on NCBI RefSeq
#' compared to other sources such as Ensembl and GENCODE, which track a stable
#' release as the latest "reference" assembly.
#'
#' @export
#' @note Updated 2021-01-14.
#'
#' @param file `character(1)`.
#'   RefSeq assembly summary file or URL.
#'
#' @seealso
#' - [File format details](ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt).
#'
#' @return Named `character`.
#'
#' @examples
#' file <- pasteURL(
#'     "ftp.ncbi.nlm.nih.gov",
#'     "genomes",
#'     "refseq",
#'     "vertebrate_mammalian",
#'     "Homo_sapiens",
#'     "assembly_summary.txt",
#'     protocol = "ftp"
#' )
#' x <- getRefSeqAssemblySummary(file)
#' names(x)

## nolint end

getRefSeqAssemblySummary <-
    function(file) {
        pattern <- "assembly_summary.txt"
        assert(
            isString(file),
            isMatchingFixed(pattern = pattern, x = basename(file))
        )
        file <- .cacheIt(file)
        lines <- import(
            file = file,
            format = "lines",
            skip = 1L,
            quiet = TRUE
        )
        names <- strsplit(
            x = sub(pattern = "^#\\s", replacement = "", x = lines[[1L]]),
            split = "\\t"
        )[[1L]]
        values <- strsplit(x = lines[[2L]], split = "\\t")[[1L]]
        x <- as.character(values[seq_len(20L)])
        names(x) <- names[seq_len(20L)]
        x <- x[nzchar(x)]
        x
    }
