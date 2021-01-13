#' Get the RefSeq assembly metadata
#'
#' @note Updated 2021-01-08.
#' @noRd
#'
#' @param baseURL `character(1)`.
#'   RefSeq organism base URL.
#'   (e.g. "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/
#'       "vertebrate_mammalian/Homo_sapiens).
#'
#' @seealso
#' - ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt
#'
#' @return Named `character`.
.getRefSeqAssemblySummary <-
    function(baseURL) {
        lines <- import(
            file = pasteURL(baseURL, "assembly_summary.txt"),
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



#' Get the RefSeq base genome URL for an organism
#'
#' @note Updated 2021-01-08.
#' @noRd
.getRefSeqGenomeURL <- function(
    organism,
    taxonomicGroup = NULL,
    quiet = FALSE
) {
    assert(
        isOrganism(organism),
        isString(taxonomicGroup, nullOK = TRUE),
        isFlag(quiet)
    )
    if (isFALSE(quiet)) {
        alert(sprintf(
            "Locating {.emph %s} genome on RefSeq FTP server.",
            organism
        ))
    }
    baseURL <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq"
    if (is.null(taxonomicGroup)) {
        taxonomicGroups <- getURLDirList(url = baseURL)
        keep <- grepl(pattern = "^[a-z_]+$", x = taxonomicGroups)
        taxonomicGroups <- sort(taxonomicGroups[keep])
        list <- bplapply(
            X = taxonomicGroups,
            baseURL = baseURL,
            FUN = function(taxonomicGroup, baseURL) {
                url <- pasteURL(baseURL, taxonomicGroup)
                x <- getURLDirList(url = url)
                keep <- grepl(pattern = "^[A-Z][a-z]+_[a-z]+$", x = x)
                x <- sort(x[keep])
                x
            }
        )
        names(list) <- taxonomicGroups
        match <- vapply(
            X = list,
            organism = gsub(pattern = " ", replacement = "_", x = organism),
            FUN = function(strings, organism) {
                isSubset(x = organism, y = strings)
            },
            FUN.VALUE = logical(1L),
            USE.NAMES = TRUE
        )
        taxonomicGroup <- names(match)[match]
        assert(isString(taxonomicGroup))
    }
    url <- pasteURL(
        baseURL,
        taxonomicGroup,
        gsub(pattern = " ", replacement = "_", x = organism)
    )
    if (isFALSE(quiet)) {
        dl(c("URL" = url))
    }
    url
}
