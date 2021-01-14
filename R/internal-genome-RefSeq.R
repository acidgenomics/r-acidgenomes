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
