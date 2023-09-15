#' Map Ensembl release to archive URL.
#'
#' @note Updated 2023-04-26.
#' @export
#'
#' @param release `integer(1)` or `character(1)`.
#' Ensembl release (e.g. `100`).
#'
#' @return `character(1)`.
#' URL.
#'
#' @seealso
#' - `biomaRt::listEnsemblArchives()`.
#'
#' @examples
#' tryCatch(
#'     expr = mapEnsemblReleaseToURL(release = 100L),
#'     error = function(e) message(e)
#' )
mapEnsemblReleaseToURL <- function(release) {
    assert(requireNamespaces("biomaRt"))
    currentURL <- "https://useast.ensembl.org"
    if (is.null(release)) {
        return(currentURL)
    }
    release <- as.character(release)
    assert(isString(release))
    map <- tryCatch(
        expr = biomaRt::listEnsemblArchives(),
        error = function(e) {
            abort(sprintf(
                "{.pkg %s}::{.fun %s} error.",
                "biomaRt", "listEnsemblArchives"
            ))
        }
    )
    assert(
        is.data.frame(map),
        isSubset(c("url", "version"), colnames(map))
    )
    if (!release %in% map[["version"]]) {
        abort(sprintf(
            "Supported Ensembl releases: %s.",
            toInlineString(map[["version"]], n = 10L)
        ))
    }
    ## Extract the matching row, so we can check if releast is current.
    which <- match(x = release, table = map[["version"]])
    x <- map[which, , drop = FALSE]
    isCurrent <- identical(x[1L, "current_release"], "*")
    if (isTRUE(isCurrent)) {
        return(currentURL)
    }
    url <- x[1L, "url"]
    assert(isTRUE(grepl("ensembl\\.org", url)))
    url
}
